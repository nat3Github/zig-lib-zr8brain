const std = @import("std");
const zpffft = @import("zpffft");

/// Cache for RealFFT(f64) setups keyed by len_bits. Thread-unsafe; FIRFilterCache.getFFT locks it.
pub const RealFFTCache = struct {
    const max_len_bits = 32;
    entries: [max_len_bits]?*RealFFT(f64) = [_]?*RealFFT(f64){null} ** max_len_bits,
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator) RealFFTCache {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *RealFFTCache) void {
        for (&self.entries) |*e| {
            if (e.*) |fft| {
                fft.deinit();
                e.* = null;
            }
        }
    }

    pub fn get(self: *RealFFTCache, len_bits: usize) !*RealFFT(f64) {
        if (len_bits >= max_len_bits) return error.LenBitsTooLarge;
        if (self.entries[len_bits]) |cached| return cached;
        const fft = try RealFFT(f64).init(self.allocator, len_bits);
        self.entries[len_bits] = fft;
        return fft;
    }
};

pub fn RealFFT(comptime T: type) type {
    return struct {
        const Self = @This();
        const FloatEngine = if (T == f32) zpffft.Float32 else zpffft.Float64;
        
        len_bits: usize,
        len: usize,
        inv_mul_const: T,
        
        setup: *FloatEngine.Setup,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator, len_bits: usize) !*Self {
            const len = @as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(len_bits));
            
            const self = try allocator.create(Self);
            errdefer allocator.destroy(self);
            
            const setup = try FloatEngine.Setup.init(allocator, @as(i32, @intCast(len)), .Real);
            errdefer setup.deinit();

            self.* = .{
                .len_bits = len_bits,
                .len = len,
                .inv_mul_const = 1.0 / @as(T, @floatFromInt(len)),
                .setup = setup,
                .allocator = allocator,
            };

            return self;
        }

        pub fn deinit(self: *Self) void {
            self.setup.deinit();
            self.allocator.destroy(self);
        }

        /// `work`: caller-owned scratch of at least len * 2 (see allocWork). The FFT itself is
        /// read-only, so one instance can serve many threads, each with its own work.
        pub fn forward(self: *const Self, p: []T, work: []align(64) T) void {
            self.run(p, work, .Forward);
        }

        pub fn inverse(self: *const Self, p: []T, work: []align(64) T) void {
            self.run(p, work, .Backward);
        }

        pub fn allocWork(self: *const Self, allocator: std.mem.Allocator) ![]align(64) T {
            return allocator.alignedAlloc(T, .@"64", self.len * 2);
        }

        fn run(self: *const Self, p: []T, work: []align(64) T, direction: zpffft.Direction) void {
            std.debug.assert(p.len >= self.len and work.len >= self.len * 2);
            const out = work[0..self.len];
            const scr = work[self.len..][0..self.len];
            // Cannot fail: out and scr are disjoint 64-aligned halves of work (alignment() <= 64),
            // len floats each; p is outside work.
            if (std.mem.isAligned(@intFromPtr(p.ptr), self.setup.alignment())) {
                // In place: no copy back (BlockConvolver's blocks are 64-aligned).
                FloatEngine.transformOrdered(self.setup, p[0..self.len], p[0..self.len], scr, direction) catch unreachable;
                return;
            }
            // An unaligned p (plain alloc) is copied to out and transformed in place there.
            @memcpy(out, p[0..self.len]);
            FloatEngine.transformOrdered(self.setup, out, out, scr, direction) catch unreachable;
            @memcpy(p[0..self.len], out);
        }

        pub fn convertToZP(self: *const Self, p: []T) void {
            _ = self;
            var i: usize = 2;
            while (i < p.len) : (i += 2) {
                p[i + 1] = p[i];
            }
        }

        pub fn multiplyBlocks(self: *const Self, ip1: []const T, op: []T) void {
            _ = self;
            std.debug.assert(ip1.len == op.len);
            
            op[0] *= ip1[0];
            op[1] *= ip1[1];

            var i: usize = 2;
            while (i < op.len) : (i += 2) {
                const t = op[i] * ip1[i] - op[i + 1] * ip1[i + 1];
                op[i + 1] = op[i] * ip1[i + 1] + op[i + 1] * ip1[i];
                op[i] = t;
            }
        }

        pub fn multiplyBlocksZP(self: *const Self, ip: []const T, op: []T) void {
            // "ip" contains zero-phase response
            const len = self.len;
            std.debug.assert(ip.len >= len and op.len >= len);
            
            const VectorSize = 4;
            const V = @Vector(VectorSize, T);
            
            var i: usize = 0;
            const limit = len - (len % VectorSize);
            while (i < limit) : (i += VectorSize) {
                const vec_ip: V = ip[i..][0..VectorSize].*;
                const vec_op: V = op[i..][0..VectorSize].*;
                const result = vec_ip * vec_op;
                @memcpy(op[i..][0..VectorSize], &@as([VectorSize]T, result));
            }
            
            while (i < len) : (i += 1) {
                op[i] *= ip[i];
            }
        }
    };
}
