const std = @import("std");
const base = @import("base.zig");
pub const halfband_data = @import("halfband_data.zig");

pub fn getHBFilter(req_atten: f64, steep_index: i32) halfband_data.HBFilter {
    const groups = [_][]const halfband_data.HBFilter{
        &halfband_data.FltGroupA,
        &halfband_data.FltGroupB,
        &halfband_data.FltGroupC,
        &halfband_data.FltGroupD,
        &halfband_data.FltGroupE,
        &halfband_data.FltGroupF,
        &halfband_data.FltGroupG,
    };
    
    var idx = @as(usize, @intCast(@max(0, steep_index)));
    if (idx >= groups.len) idx = groups.len - 1;
    const group = groups[idx];
    
    var k: usize = 0;
    while (k != group.len - 1 and group[k].atten < req_atten) {
        k += 1;
    }
    return group[k];
}

pub fn getHBFilterThird(req_atten: f64, steep_index: i32) halfband_data.HBFilter {
    const groups = [_][]const halfband_data.HBFilter{
        &halfband_data.FltGroupA_Third,
        &halfband_data.FltGroupB_Third,
        &halfband_data.FltGroupC_Third,
        &halfband_data.FltGroupD_Third,
        &halfband_data.FltGroupE_Third,
        &halfband_data.FltGroupF_Third,
        &halfband_data.FltGroupG_Third,
    };
    
    var idx = @as(usize, @intCast(@max(0, steep_index)));
    if (idx >= groups.len) idx = groups.len - 1;
    const group = groups[idx];
    
    var k: usize = 0;
    while (k != group.len - 1 and group[k].atten < req_atten) {
        k += 1;
    }
    return group[k];
}

pub fn HBUpsampler(comptime T: type) type {
    return struct {
        const Self = @This();

        const buf_len_bits = 9;
        const buf_len = 1 << buf_len_bits;
        const buf_len_mask = buf_len - 1;

        buf: [buf_len + 27]T = undefined,
        fltp: []const f64,
        fll: usize,
        fl2: usize,
        flo: usize,
        flb: usize,
        latency_frac: T,
        latency: usize,
        latency_left: usize,
        buf_left: usize,
        write_pos: usize,
        read_pos: usize,
        do_consume_latency: bool,
        buf_rp_offset: usize,
        is_third: bool,
        steep_index: i32,

        pub fn init(req_atten: f64, steep_index: i32, is_third: bool, prev_latency: T, do_consume_latency: bool) Self {
            const filter = if (is_third) getHBFilterThird(req_atten, steep_index) else getHBFilter(req_atten, steep_index);
            
            const fltt = filter.taps.len;
            var lat_frac = prev_latency * 2.0;
            const lat = @as(usize, @intFromFloat(lat_frac));
            lat_frac -= @as(T, @floatFromInt(lat));

            var self = Self{
                .fltp = filter.taps,
                .fll = fltt - 1,
                .fl2 = fltt,
                .flo = (fltt - 1) + fltt,
                .flb = 0,
                .latency_frac = lat_frac,
                .latency = lat,
                .latency_left = 0,
                .buf_left = 0,
                .write_pos = 0,
                .read_pos = 0,
                .do_consume_latency = do_consume_latency,
                .buf_rp_offset = fltt - 1,
                .is_third = is_third,
                .steep_index = steep_index,
            };

            if (do_consume_latency) {
                self.flb = buf_len - self.fll;
            } else {
                self.latency += self.fl2 + self.fl2;
                self.flb = buf_len - self.flo;
            }

            self.clear();
            return self;
        }

        pub fn getInLenBeforeOutPos(self: Self, req_out_pos: usize) usize {
            return self.fl2 + @as(usize, @intFromFloat((@as(T, @floatFromInt(self.latency)) + self.latency_frac + @as(T, @floatFromInt(req_out_pos))) * 0.5));
        }

        pub fn getLatency(self: Self) usize {
            return self.latency;
        }

        pub fn getLatencyFrac(self: *const Self) f64 {
            return @as(f64, @floatCast(self.latency_frac));
        }

        pub fn getName(ptr: *anyopaque) []const u8 {
            const self: *Self = @ptrCast(@alignCast(ptr));
            return if (self.is_third) "HBUpsampler(3x)" else "HBUpsampler(2x)";
        }

        pub fn getMaxOutLen(self: *const Self, max_in_len: usize) usize {
            _ = self;
            return max_in_len * 2;
        }

        pub fn processor(self: *Self) base.Processor {
            return .{
                .ptr = self,
                .vtable = &.{
                    .getName = @ptrCast(&getName),
                    .getInLenBeforeOutPos = @ptrCast(&getInLenBeforeOutPos),
                    .getLatency = @ptrCast(&getLatency),
                    .getLatencyFrac = @ptrCast(&getLatencyFrac),
                    .getMaxOutLen = @ptrCast(&getMaxOutLen),
                    .clear = @ptrCast(&clear),
                    .process = @ptrCast(&processGeneric),
                    .deinit = @ptrCast(&deinitGeneric),
                },
            };
        }

        fn processGeneric(ptr: *anyopaque, ip: [*]const T, l: usize, op0: *[*]T) usize {
            const self: *Self = @ptrCast(@alignCast(ptr));
            const written = self.process(ip[0..l], op0.*[0 .. l * 2]);
            return written;
        }

        fn deinitGeneric(ptr: *anyopaque, allocator: std.mem.Allocator) void {
            const self: *Self = @ptrCast(@alignCast(ptr));
            allocator.destroy(self);
        }

        pub fn clear(self: *Self) void {
            if (self.do_consume_latency) {
                self.latency_left = self.latency;
                self.buf_left = 0;
            } else {
                self.latency_left = 0;
                self.buf_left = self.fl2;
            }

            self.write_pos = 0;
            self.read_pos = self.flb;
            @memset(&self.buf, 0);
        }

        inline fn convolve(self: *Self, op: *[]T, limit: usize) void {
            var rpos = self.read_pos;
            var op_idx: usize = 0;
            while (op_idx < limit) : (op_idx += 2) {
                const rp_idx = self.buf_rp_offset + rpos;
                op.*[op_idx] = self.buf[rp_idx];
                op.*[op_idx + 1] = base.firSymDot(T, self.fltp, &self.buf, rp_idx);

                rpos = (rpos + 1) & buf_len_mask;
            }
            self.read_pos = rpos;
        }

        pub fn process(self: *Self, ip_ptr: []const T, op_full: []T) usize {
            var l = ip_ptr.len;
            var ip = ip_ptr;
            var op = op_full;

            while (l > 0) {
                const b = @min(l, @min(buf_len - self.write_pos, self.flb - self.buf_left));

                @memcpy(self.buf[self.write_pos .. self.write_pos + b], ip[0 .. b]);

                const ec = @as(isize, @intCast(self.flo)) - @as(isize, @intCast(self.write_pos));
                if (ec > 0) {
                    const copy_len = @min(b, @as(usize, @intCast(ec)));
                    @memcpy(self.buf[self.write_pos + buf_len .. self.write_pos + buf_len + copy_len], ip[0 .. copy_len]);
                }

                ip = ip[b..];
                self.write_pos = (self.write_pos + b) & buf_len_mask;
                l -= b;
                self.buf_left += b;

                // Produce output
                const c = @as(isize, @intCast(self.buf_left)) - @as(isize, @intCast(self.fl2));
                if (c > 0) {
                    const items = @as(usize, @intCast(c)) * 2;
                    self.convolve(&op, items);
                    op = op[items..];
                    self.buf_left -= @as(usize, @intCast(c));
                }
            }

            var total_out = op_full.len - op.len;

            // Consume latency from front of output (only active when do_consume_latency=true)
            if (self.latency_left > 0) {
                if (self.latency_left >= total_out) {
                    self.latency_left -= total_out;
                    return 0;
                }
                total_out -= self.latency_left;
                std.mem.copyForwards(T, op_full[0..total_out], op_full[self.latency_left .. self.latency_left + total_out]);
                self.latency_left = 0;
            }

            return total_out;
        }
    };
}

pub fn HBDownsampler(comptime T: type) type {
    return struct {
        const Self = @This();

        const buf_len_bits = 10;
        const buf_len = 1 << buf_len_bits;
        const buf_len_mask = buf_len - 1;

        buf1: [buf_len + 27]T = undefined,
        buf2: [buf_len + 27]T = undefined,
        fltp: []const f64,
        fll: usize,
        fl2: usize,
        flo: usize,
        flb: usize,
        latency_frac: T,
        latency: usize,
        latency_left: usize,
        buf_left: usize,
        write_pos1: usize,
        write_pos2: usize,
        read_pos: usize,
        buf_rp1_offset: usize,
        buf_rp2_offset: usize,
        is_third: bool,
        steep_index: i32,

        pub fn init(req_atten: f64, steep_index: i32, is_third: bool, prev_latency: T) Self {
            const filter = if (is_third) getHBFilterThird(req_atten, steep_index) else getHBFilter(req_atten, steep_index);
            
            const fltt = filter.taps.len;
            var lat_frac = prev_latency * 0.5;
            const lat = @as(usize, @intFromFloat(lat_frac));
            lat_frac -= @as(T, @floatFromInt(lat));

            var self = Self{
                .fltp = filter.taps,
                .fll = fltt,
                .fl2 = fltt - 1,
                .flo = fltt + (fltt - 1),
                .flb = buf_len - fltt,
                .latency_frac = lat_frac,
                .latency = lat,
                .latency_left = 0,
                .buf_left = 0,
                .write_pos1 = 0,
                .write_pos2 = 0,
                .read_pos = buf_len - fltt,
                .buf_rp1_offset = fltt,
                .buf_rp2_offset = fltt - 1,
                .is_third = is_third,
                .steep_index = steep_index,
            };

            self.clear();
            return self;
        }

        pub fn getInLenBeforeOutPos(self: Self, req_out_pos: usize) usize {
            return self.flo + @as(usize, @intFromFloat((@as(T, @floatFromInt(self.latency)) + self.latency_frac + @as(T, @floatFromInt(req_out_pos))) * 2.0));
        }

        pub fn getLatency(self: Self) usize {
            return self.latency;
        }

        pub fn getLatencyFrac(self: *const Self) f64 {
            return @as(f64, @floatCast(self.latency_frac));
        }

        pub fn getName(ptr: *anyopaque) []const u8 {
            const self: *Self = @ptrCast(@alignCast(ptr));
            return if (self.is_third) "HBDownsampler(3x)" else "HBDownsampler(2x)";
        }

        pub fn getMaxOutLen(self: *const Self, max_in_len: usize) usize {
            _ = self;
            return (max_in_len + 1) / 2;
        }

        pub fn processor(self: *Self) base.Processor {
            return .{
                .ptr = self,
                .vtable = &.{
                    .getName = @ptrCast(&getName),
                    .getInLenBeforeOutPos = @ptrCast(&getInLenBeforeOutPos),
                    .getLatency = @ptrCast(&getLatency),
                    .getLatencyFrac = @ptrCast(&getLatencyFrac),
                    .getMaxOutLen = @ptrCast(&getMaxOutLen),
                    .clear = @ptrCast(&clear),
                    .process = @ptrCast(&processGeneric),
                    .deinit = @ptrCast(&deinitGeneric),
                },
            };
        }

        fn processGeneric(ptr: *anyopaque, ip: [*]const T, l: usize, op0: *[*]T) usize {
            const self: *Self = @ptrCast(@alignCast(ptr));
            const written = self.process(ip[0..l], op0.*[0 .. (l + 1) / 2]);
            return written;
        }

        fn deinitGeneric(ptr: *anyopaque, allocator: std.mem.Allocator) void {
            const self: *Self = @ptrCast(@alignCast(ptr));
            allocator.destroy(self);
        }

        pub fn clear(self: *Self) void {
            self.latency_left = self.latency;
            self.buf_left = 0;
            self.write_pos1 = 0;
            self.write_pos2 = 0;
            self.read_pos = self.flb;
            
            @memset(&self.buf1, 0);
            @memset(&self.buf2, 0);
        }

        inline fn convolve(self: *Self, op: *[]T, limit: usize) void {
            var rpos = self.read_pos;
            var op_idx: usize = 0;
            while (op_idx < limit) : (op_idx += 1) {
                const rp1_idx = self.buf_rp1_offset + rpos;
                const rp2_idx = self.buf_rp2_offset + rpos;

                op.*[op_idx] = self.buf1[rp1_idx] + base.firSymDot(T, self.fltp, &self.buf2, rp2_idx);

                rpos = (rpos + 1) & buf_len_mask;
            }
            self.read_pos = rpos;
        }

        pub fn process(self: *Self, ip_ptr: []const T, op_full: []T) usize {
            var l = ip_ptr.len;
            var ip = ip_ptr;
            var op = op_full;

            while (l > 0) {
                if (self.write_pos1 != self.write_pos2) {
                    self.buf2[self.write_pos2] = ip[0];
                    if (self.write_pos2 < self.flo) {
                        self.buf2[self.write_pos2 + buf_len] = ip[0];
                    }
                    ip = ip[1..];
                    self.write_pos2 = self.write_pos1;
                    l -= 1;
                    self.buf_left += 1;
                }

                const b1 = @min((l + 1) >> 1, @min(buf_len - self.write_pos1, self.flb - self.buf_left));
                const b2 = if (b1 * 2 > l) b1 - 1 else b1;

                var ip_idx: usize = 0;
                var wp_idx: usize = 0;
                while (ip_idx < b2 * 2) {
                    self.buf1[self.write_pos1 + wp_idx] = ip[ip_idx];
                    self.buf2[self.write_pos1 + wp_idx] = ip[ip_idx + 1];
                    wp_idx += 1;
                    ip_idx += 2;
                }

                if (b1 != b2) {
                    self.buf1[self.write_pos1 + wp_idx] = ip[ip_idx];
                    ip_idx += 1;
                }
                ip = ip[ip_idx..];

                const ec = @as(isize, @intCast(self.flo)) - @as(isize, @intCast(self.write_pos1));
                if (ec > 0) {
                    const ec_u = @as(usize, @intCast(ec));
                    @memcpy(self.buf1[self.write_pos1 + buf_len .. self.write_pos1 + buf_len + @min(b1, ec_u)], self.buf1[self.write_pos1 .. self.write_pos1 + @min(b1, ec_u)]);
                    @memcpy(self.buf2[self.write_pos1 + buf_len .. self.write_pos1 + buf_len + @min(b2, ec_u)], self.buf2[self.write_pos1 .. self.write_pos1 + @min(b2, ec_u)]);
                }

                self.write_pos1 = (self.write_pos1 + b1) & buf_len_mask;
                self.write_pos2 = (self.write_pos2 + b2) & buf_len_mask;
                l -= b1 + b2;
                self.buf_left += b2;

                // Produce output
                const c = @as(isize, @intCast(self.buf_left)) - @as(isize, @intCast(self.fl2));
                if (c > 0) {
                    const items = @as(usize, @intCast(c));
                    self.convolve(&op, items);
                    op = op[items..];
                    self.buf_left -= items;
                }
            }

            const ol = op_full.len - op.len;
            return ol;
        }
    };
}
