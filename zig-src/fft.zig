const std = @import("std");
const zpffft = @import("zpffft");

/// Cache for RealFFT(f64) setups keyed by len_bits. Thread-unsafe; caller owns.
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
        work: []align(64) T,
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
                .work = try allocator.alignedAlloc(T, .@"64", len * 2),
                .allocator = allocator,
            };
            @memset(self.work, 0);
            
            return self;
        }

        pub fn deinit(self: *Self) void {
            self.setup.deinit();
            self.allocator.free(self.work);
            self.allocator.destroy(self);
        }

        pub fn forward(self: *const Self, p: []T) void {
            std.debug.assert(p.len >= self.len);
            const out = self.work[0..self.len];
            const scr = self.work[self.len..];
            FloatEngine.transformOrdered(self.setup, p[0..self.len], out, scr, .Forward);
            @memcpy(p[0..self.len], out);
        }

        pub fn inverse(self: *const Self, p: []T) void {
            std.debug.assert(p.len >= self.len);
            const out = self.work[0..self.len];
            const scr = self.work[self.len..];
            FloatEngine.transformOrdered(self.setup, p[0..self.len], out, scr, .Backward);
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
