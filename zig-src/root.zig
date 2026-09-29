const std = @import("std");
pub const base = @import("base.zig");
pub const fft = @import("fft.zig");
pub const filters = @import("filters.zig");
pub const halfband = @import("halfband.zig");
pub const interpolator = @import("interpolator.zig");
pub const resampler = @import("resampler.zig");

pub const FIRFilterCache = filters.FIRFilterCache;
pub const RealFFTCache = @import("fft.zig").RealFFTCache;
pub const FracDelayFilterBankCache = @import("interpolator.zig").FracDelayFilterBankCache;

/// Quality presets matching the C++ CDSPResampler named variants.
pub const ResamplerQuality = enum {
    /// CDSPResampler16: 16-bit quality (~136 dB attenuation)
    quality16,
    /// CDSPResampler16IR: 16-bit quality, intermediate rate (~110 dB attenuation)
    quality16ir,
    /// CDSPResampler24: 24-bit quality (~180 dB attenuation) — default
    quality24,

    pub fn attenuation(self: ResamplerQuality) f64 {
        return switch (self) {
            .quality16 => 136.45,
            .quality16ir => 109.56,
            .quality24 => 180.15,
        };
    }
};

pub const R8bResampler = struct {
    inner: *resampler.Resampler(f64),
    allocator: std.mem.Allocator,
    input_conv: []f64,

    pub fn init(allocator: std.mem.Allocator, input_rate: f64, output_rate: f64, max_in_len: u32, req_trans_band: f64, cache: *FIRFilterCache) !R8bResampler {
        return initWithQuality(allocator, input_rate, output_rate, max_in_len, req_trans_band, .quality24, cache);
    }

    pub fn initWithQuality(allocator: std.mem.Allocator, input_rate: f64, output_rate: f64, max_in_len: u32, req_trans_band: f64, quality: ResamplerQuality, cache: *FIRFilterCache) !R8bResampler {
        const inner = try resampler.Resampler(f64).init(allocator, input_rate, output_rate, max_in_len, req_trans_band, quality.attenuation(), .linearPhase, cache);
        errdefer inner.deinit();

        const input_conv = try allocator.alloc(f64, max_in_len);

        return .{
            .inner = inner,
            .allocator = allocator,
            .input_conv = input_conv,
        };
    }

    pub fn deinit(self: *R8bResampler) void {
        self.inner.deinit();
        self.allocator.free(self.input_conv);
    }

    pub fn clear(self: *R8bResampler) void {
        self.inner.clear();
    }

    pub fn process(self: *R8bResampler, input: []const f32, output: []f32) usize {
        const l = @min(input.len, self.input_conv.len);
        convert(self.input_conv[0..l], input[0..l]);
        const out = self.inner.processRaw(self.input_conv[0..l]);
        const out_len = @min(out.len, output.len);
        convert(output[0..out_len], out[0..out_len]);
        return out_len;
    }

    /// Explicit SIMD: LLVM leaves these f32<->f64 loops scalar.
    fn convert(dst: anytype, src: anytype) void {
        const D = @typeInfo(@TypeOf(dst)).pointer.child;
        const S = @typeInfo(@TypeOf(src)).pointer.child;
        const n = 8;
        var i: usize = 0;
        while (i + n <= dst.len) : (i += n) dst[i..][0..n].* = @as(@Vector(n, D), @floatCast(@as(@Vector(n, S), src[i..][0..n].*)));
        while (i < dst.len) : (i += 1) dst[i] = @floatCast(src[i]);
    }

    pub fn getInLenBeforeOutPos(self: *const R8bResampler, req_out_pos: i32) i32 {
        return @intCast(self.inner.getInLenBeforeOutPos(@intCast(req_out_pos)));
    }

    pub fn getLatency(self: *const R8bResampler) i32 {
        return @intCast(self.inner.getLatency());
    }

    pub fn getLatencyFrac(self: *const R8bResampler) f64 {
        return self.inner.latency_frac;
    }

    pub fn getMaxOutLen(self: *const R8bResampler, max_in_len: i32) i32 {
        return @intCast(self.inner.getMaxOutLen(@intCast(max_in_len)));
    }
};
