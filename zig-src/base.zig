const std = @import("std");

/// Lane count for the FIR inner loops. Float adds are not reassociated by the
/// compiler, so a scalar `s += f*x` loop is one serial add chain; vector
/// partial sums run `lanes` chains in parallel (C++ r8brain uses SSE2/NEON here).
pub const fir_lanes = 4;

/// Σ flt[i] * x[i] for i < flt.len.
pub inline fn firDot(comptime T: type, flt: []const f64, x: []const T) T {
    const V = @Vector(fir_lanes, T);
    var acc: V = @splat(0);
    var i: usize = 0;
    while (i + fir_lanes <= flt.len) : (i += fir_lanes) {
        const f: V = @floatCast(@as(@Vector(fir_lanes, f64), flt[i..][0..fir_lanes].*));
        acc += f * @as(V, x[i..][0..fir_lanes].*);
    }
    var s = @reduce(.Add, acc);
    while (i < flt.len) : (i += 1) s += @as(T, @floatCast(flt[i])) * x[i];
    return s;
}

/// Symmetric half-band FIR: Σ flt[i] * (x[c + 1 + i] + x[c - i]). Needs c + 1 >= flt.len.
pub inline fn firSymDot(comptime T: type, flt: []const f64, x: []const T, c: usize) T {
    const V = @Vector(fir_lanes, T);
    const rev = comptime blk: {
        var r: [fir_lanes]i32 = undefined;
        for (&r, 0..) |*e, k| e.* = fir_lanes - 1 - @as(i32, @intCast(k));
        break :blk r;
    };
    var acc: V = @splat(0);
    var i: usize = 0;
    while (i + fir_lanes <= flt.len) : (i += fir_lanes) {
        const f: V = @floatCast(@as(@Vector(fir_lanes, f64), flt[i..][0..fir_lanes].*));
        const fwd: V = x[c + 1 + i ..][0..fir_lanes].*;
        const back_raw: V = x[c + 1 - i - fir_lanes ..][0..fir_lanes].*;
        acc += f * (fwd + @shuffle(T, back_raw, undefined, rev));
    }
    var s = @reduce(.Add, acc);
    while (i < flt.len) : (i += 1) s += @as(T, @floatCast(flt[i])) * (x[c + 1 + i] + x[c - i]);
    return s;
}

/// Sine signal generator without biasing.
pub const SineGen = struct {
    svalue1: f64,
    svalue2: f64,
    sincr: f64,

    pub fn init(si: f64, ph: f64, gain: f64) SineGen {
        return .{
            .svalue1 = @sin(ph) * gain,
            .svalue2 = @sin(ph - si) * gain,
            .sincr = 2.0 * @cos(si),
        };
    }

    pub fn generate(self: *SineGen) f64 {
        const res = self.svalue1;
        self.svalue1 = self.sincr * res - self.svalue2;
        self.svalue2 = res;
        return res;
    }
};

/// Base interface for DSP processors.
pub const Processor = struct {
    ptr: *anyopaque,
    vtable: *const VTable,

    pub const VTable = struct {
        getName: *const fn (ptr: *anyopaque) []const u8,
        getInLenBeforeOutPos: *const fn (ptr: *anyopaque, req_out_pos: usize) usize,
        getLatency: *const fn (ptr: *anyopaque) usize,
        getLatencyFrac: *const fn (ptr: *anyopaque) f64,
        getMaxOutLen: *const fn (ptr: *anyopaque, max_in_len: usize) usize,
        clear: *const fn (ptr: *anyopaque) void,
        process: *const fn (ptr: *anyopaque, ip: [*]const f64, l: usize, op0: *[*]f64) usize,
        deinit: *const fn (ptr: *anyopaque, allocator: std.mem.Allocator) void,
    };

    pub fn getName(self: Processor) []const u8 {
        return self.vtable.getName(self.ptr);
    }

    pub fn getInLenBeforeOutPos(self: Processor, req_out_pos: usize) usize {
        return self.vtable.getInLenBeforeOutPos(self.ptr, req_out_pos);
    }
    pub fn getLatency(self: Processor) usize {
        return self.vtable.getLatency(self.ptr);
    }
    pub fn getLatencyFrac(self: Processor) f64 {
        return self.vtable.getLatencyFrac(self.ptr);
    }
    pub fn getMaxOutLen(self: Processor, max_in_len: usize) usize {
        return self.vtable.getMaxOutLen(self.ptr, max_in_len);
    }
    pub fn clear(self: Processor) void {
        self.vtable.clear(self.ptr);
    }
    pub fn process(self: Processor, ip: [*]const f64, l: usize, op0: *[*]f64) usize {
        return self.vtable.process(self.ptr, ip, l, op0);
    }
    pub fn deinit(self: Processor, allocator: std.mem.Allocator) void {
        self.vtable.deinit(self.ptr, allocator);
    }
};

/// Helper function to calculate bit occupancy (equivalent to 32 - clz).
pub fn getBitOccupancy(v: usize) u8 {
    if (v == 0) return 0;
    // @clz works for u32 or u64, we'll use u64 if needed, but the original used 32
    return @as(u8, @intCast(32 - @clz(@as(u32, @intCast(v)))));
}

pub fn normalizeFIRFilter(p: []f64, l: usize, dc_gain: f64, pstep: usize) void {
    var s: f64 = 0;
    var i: usize = 0;
    while (i < l) : (i += 1) {
        s += p[i * pstep];
    }
    const mul = if (s == 0) 0 else dc_gain / s;
    i = 0;
    while (i < l) : (i += 1) {
        p[i * pstep] *= mul;
    }
}

pub fn calcSpline2p8Coeffs(c: []f64, xm3: f64, xm2: f64, xm1: f64, x0: f64, x1: f64, x2: f64, x3: f64, x4: f64) void {
    c[0] = x0;
    c[1] = (61.0 * (x1 - xm1) + 16.0 * (xm2 - x2) + 3.0 * (x3 - xm3)) * 1.31578947368421052e-2;
    c[2] = (106.0 * (xm1 + x1) + 10.0 * x3 + 6.0 * xm3 - 3.0 * x4 - 29.0 * (xm2 + x2) - 167.0 * x0) * 1.31578947368421052e-2;
}

pub fn calcSpline3p8Coeffs(c: []f64, xm3: f64, xm2: f64, xm1: f64, x0: f64, x1: f64, x2: f64, x3: f64, x4: f64) void {
    c[0] = x0;
    c[1] = (61.0 * (x1 - xm1) + 16.0 * (xm2 - x2) + 3.0 * (x3 - xm3)) * 1.31578947368421052e-2;
    c[2] = (106.0 * (xm1 + x1) + 10.0 * x3 + 6.0 * xm3 - 3.0 * x4 - 29.0 * (xm2 + x2) - 167.0 * x0) * 1.31578947368421052e-2;
    c[3] = (91.0 * (x0 - x1) + 45.0 * (x2 - xm1) + 13.0 * (xm2 - x3) + 3.0 * (x4 - xm3)) * 1.31578947368421052e-2;
}
