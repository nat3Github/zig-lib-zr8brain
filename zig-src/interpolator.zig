const std = @import("std");
const base = @import("base.zig");
const filters = @import("filters.zig");

pub const FracDelayFilterBank = struct {
    const Self = @This();

    filter_len: usize,
    filter_fracs: usize,
    element_size: usize,
    interp_points: usize,
    req_atten: f64,
    is_third: bool,
    filter_size: usize,
    table: []align(64) f64,
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, filter_fracs: i32, element_size: usize, interp_points: usize, req_atten_in: f64, is_third: bool, correction_gain: f64) !*Self {
        var req_atten = req_atten_in;
        var flt_len: i32 = 0;
        const params = getWinParams(&req_atten, is_third, &flt_len);

        const self = try allocator.create(Self);
        errdefer allocator.destroy(self);

        const actual_fracs = if (filter_fracs == -1) 
            @as(usize, @intFromFloat(@ceil(std.math.pow(f64, 6.4, req_atten / 50.0)))) 
            else @as(usize, @intCast(filter_fracs));

        const filter_size = @as(usize, @intCast(flt_len)) * element_size;
        const table_len = filter_size * (actual_fracs + interp_points);
        const table = try allocator.alignedAlloc(f64, .@"64", table_len);
        @memset(table, 0);

        self.* = .{
            .filter_len = @as(usize, @intCast(flt_len)),
            .filter_fracs = actual_fracs,
            .element_size = element_size,
            .interp_points = interp_points,
            .req_atten = req_atten,
            .is_third = is_third,
            .filter_size = filter_size,
            .table = table,
            .allocator = allocator,
        };

        var sinc = filters.SincFilterGen(f64){ .len2 = @as(f64, @floatFromInt(flt_len)) / 2.0, .freq2 = std.math.pi };
        const pc2 = interp_points / 2;
        
        var p_idx: usize = 0;
        var i: isize = -@as(isize, @intCast(pc2)) + 1;
        while (i <= @as(isize, @intCast(actual_fracs)) + @as(isize, @intCast(pc2))) : (i += 1) {
            sinc.frac_delay = @as(f64, @floatFromInt(@as(isize, @intCast(actual_fracs)) - i)) / @as(f64, @floatFromInt(actual_fracs));
            sinc.initFrac(.kaiser, params, true, @as(usize, @intCast(flt_len)));
            sinc.generateFrac(table[p_idx..], element_size);
            base.normalizeFIRFilter(table[p_idx..], self.filter_len, correction_gain, element_size);
            p_idx += filter_size;
        }

        const table_pos2 = filter_size;
        const table_pos3 = filter_size * 2;
        const table_pos4 = filter_size * 3;
        const table_pos5 = filter_size * 4;
        const table_pos6 = filter_size * 5;
        const table_pos7 = filter_size * 6;
        const table_pos8 = filter_size * 7;
        const table_end_idx = (actual_fracs + 1) * filter_size;
        
        var p_off: usize = 0;
        if (interp_points == 8) {
            if (element_size == 3) {
                while (p_off < table_end_idx) : (p_off += element_size) {
                    base.calcSpline2p8Coeffs(table[p_off..], table[p_off], table[p_off + table_pos2],
                        table[p_off + table_pos3], table[p_off + table_pos4], table[p_off + table_pos5],
                        table[p_off + table_pos6], table[p_off + table_pos7], table[p_off + table_pos8]);
                }
                // Like C++ shuffle2_3 under R8B_SIMD_ISH, but for fir_lanes:
                // each full block of L taps becomes [f0 x L][f1 x L][f2 x L] so
                // convolveQuadratic does contiguous loads. Tail stays interleaved.
                const L = base.fir_lanes;
                var tmp: [3 * L]f64 = undefined;
                var f: usize = 0;
                while (f < table_end_idx) : (f += filter_size) {
                    var t: usize = 0;
                    while (t + L <= self.filter_len) : (t += L) {
                        const blk = table[f + t * 3 ..][0 .. 3 * L];
                        for (0..3) |k| for (0..L) |j| {
                            tmp[k * L + j] = blk[j * 3 + k];
                        };
                        @memcpy(blk, &tmp);
                    }
                }
            } else if (element_size == 4) {
                while (p_off < table_end_idx) : (p_off += element_size) {
                    base.calcSpline3p8Coeffs(table[p_off..], table[p_off], table[p_off + table_pos2],
                        table[p_off + table_pos3], table[p_off + table_pos4], table[p_off + table_pos5],
                        table[p_off + table_pos6], table[p_off + table_pos7], table[p_off + table_pos8]);
                }
            }
        } else {
            if (element_size == 2) {
                while (p_off < table_end_idx) : (p_off += element_size) {
                    table[p_off + 1] = table[p_off + table_pos2] - table[p_off];
                }
            }
        }

        return self;
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.table);
        self.allocator.destroy(self);
    }

    pub fn getFilter(self: Self, i: usize) []const f64 {
        const start = i * self.filter_size;
        return self.table[start .. start + self.filter_size];
    }

    pub fn getWinParams(att: *f64, is_third: bool, fltlen: *i32) []const f64 {
        const Coeffs2 = [_][3]f64{
            .{ 4.1308468534586913, 1.1752580009977263, 55.5446 },
            .{ 4.4241520324148826, 1.8004881791443044, 81.4191 },
            .{ 5.2615232289173663, 1.8133318236025469, 96.3392 },
            .{ 5.9433751227216174, 1.8730186391986436, 111.1315 },
            .{ 6.8308658290513815, 1.8549555110340281, 125.4653 },
            .{ 7.6648458290312904, 1.8565766090828464, 139.7379 },
            .{ 8.2038728664307605, 1.9269521820570166, 154.0532 },
            .{ 8.7865150946655142, 1.9775307667441668, 168.2101 },
            .{ 9.5945017884101773, 1.9718456992078597, 182.1076 },
            .{ 10.5163141145985240, 1.9504067820201083, 195.5668 },
            .{ 10.2382465206362470, 2.1608923446870087, 209.0610 },
            .{ 10.9976060250714000, 2.1536533525688935, 222.5010 },
        };
        const coeffs2_base = 8;
        
        const Coeffs3 = [_][3]f64{
            .{ 3.9888564562781847, 1.5869927184268915, 66.5701 },
            .{ 4.6986694038145007, 1.8086068597928262, 86.4715 },
            .{ 5.5995071329337822, 1.8930163360942349, 106.1195 },
            .{ 6.3627287800257228, 1.9945748322093975, 125.2307 },
            .{ 7.4299550711428308, 1.9893400572347544, 144.3469 },
            .{ 8.0667715944075642, 2.0928201458699909, 163.4099 },
            .{ 8.7469970226288822, 2.1640279784268355, 181.0694 },
            .{ 10.0823430069835230, 2.0896678025321922, 199.2880 },
            .{ 10.9222206090489510, 2.1221681162186004, 216.6865 },
            .{ 21.2017743894772010, 1.1856768080118900, 233.9188 },
        };
        const coeffs3_base = 6;

        var i: usize = 0;
        if (is_third) {
            while (i != Coeffs3.len - 1 and Coeffs3[i][2] < att.*) : (i += 1) {}
            att.* = Coeffs3[i][2];
            fltlen.* = @intCast(coeffs3_base + i * 2);
            return Coeffs3[i][0..2];
        } else {
            while (i != Coeffs2.len - 1 and Coeffs2[i][2] < att.*) : (i += 1) {}
            att.* = Coeffs2[i][2];
            fltlen.* = @intCast(coeffs2_base + i * 2);
            return Coeffs2[i][0..2];
        }
    }
};

/// LRU cache for FracDelayFilterBank instances. Thread-unsafe; caller owns.
pub const FracDelayFilterBankCache = struct {
    const Self = @This();

    const Entry = struct {
        bank: *FracDelayFilterBank,
        ref_count: usize,
        filter_fracs: i32,
        element_size: usize,
        interp_points: usize,
        req_atten: f64,
        is_third: bool,
        next: ?*Entry,
    };

    allocator: std.mem.Allocator,
    head: ?*Entry,

    pub fn init(allocator: std.mem.Allocator) Self {
        return .{ .allocator = allocator, .head = null };
    }

    pub fn deinit(self: *Self) void {
        var cur = self.head;
        while (cur) |e| {
            const nxt = e.next;
            e.bank.deinit();
            self.allocator.destroy(e);
            cur = nxt;
        }
        self.head = null;
    }

    pub fn get(self: *Self, filter_fracs: i32, element_size: usize, interp_points: usize, req_atten_in: f64, is_third: bool, correction_gain: f64) !*FracDelayFilterBank {
        var rounded = req_atten_in;
        var dummy: i32 = 0;
        _ = FracDelayFilterBank.getWinParams(&rounded, is_third, &dummy);

        var prev: ?*Entry = null;
        var cur = self.head;
        while (cur) |e| {
            if (e.filter_fracs == filter_fracs and e.element_size == element_size and
                e.interp_points == interp_points and e.req_atten == rounded and e.is_third == is_third)
            {
                e.ref_count += 1;
                if (prev) |p| {
                    p.next = e.next;
                    e.next = self.head;
                    self.head = e;
                }
                return e.bank;
            }
            prev = e;
            cur = e.next;
        }

        const bank = try FracDelayFilterBank.init(self.allocator, filter_fracs, element_size, interp_points, req_atten_in, is_third, correction_gain);
        errdefer bank.deinit();

        const entry = try self.allocator.create(Entry);
        entry.* = .{ .bank = bank, .ref_count = 1, .filter_fracs = filter_fracs, .element_size = element_size, .interp_points = interp_points, .req_atten = rounded, .is_third = is_third, .next = self.head };
        self.head = entry;
        return bank;
    }

    pub fn unref(self: *Self, bank: *FracDelayFilterBank) void {
        var prev: ?*Entry = null;
        var cur = self.head;
        while (cur) |e| {
            if (e.bank == bank) {
                e.ref_count -= 1;
                if (e.ref_count == 0) {
                    if (prev) |p| p.next = e.next else self.head = e.next;
                    e.bank.deinit();
                    self.allocator.destroy(e);
                }
                return;
            }
            prev = e;
            cur = e.next;
        }
    }
};

pub fn FracInterpolator(comptime T: type) type {
    return struct {
        const Self = @This();
        const buf_len_bits = 8;
        const buf_len = 1 << buf_len_bits;
        const buf_len_mask = buf_len - 1;

        buf: [buf_len + 29]T = undefined,
        src_sample_rate: f64,
        dst_sample_rate: f64,
        init_frac_pos: f64,
        init_frac_pos_w: i32 = 0,
        latency: usize,
        latency_frac: f64,
        filter_len: usize,
        fll: usize,
        fl2: usize,
        flo: usize,
        flb: usize,
        in_step: i32 = 0,
        out_step: i32 = 0,
        latency_left: usize = 0,
        buf_left: usize = 0,
        write_pos: usize = 0,
        read_pos: usize = 0,
        in_pos_frac_w: i32 = 0,
        in_pos_frac: f64 = 0.0,
        
        frac_step: f64 = 0.0,
        in_counter: i32 = 0,
        in_pos_int: i32 = 0,
        in_pos_shift: f64 = 0.0,

        filter_bank: *FracDelayFilterBank,
        bank_cache: ?*FracDelayFilterBankCache,
        is_whole: bool,
        is_bank_owned: bool,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator, src_rate: f64, dst_rate: f64, req_atten: f64, is_third: bool, prev_latency: f64, correction_gain: f64) !*Self {
            return initWithCache(allocator, src_rate, dst_rate, req_atten, is_third, prev_latency, correction_gain, null);
        }

        pub fn initWithCache(allocator: std.mem.Allocator, src_rate: f64, dst_rate: f64, req_atten: f64, is_third: bool, prev_latency: f64, correction_gain: f64, bank_cache: ?*FracDelayFilterBankCache) !*Self {
            const self = try allocator.create(Self);
            errdefer allocator.destroy(self);

            // prev_latency is already in this stage's input samples (C++:
            // InitFracPos = PrevLatency) — no rate scaling.
            const lat = @as(usize, @intFromFloat(@floor(prev_latency)));
            const init_frac = prev_latency - @as(f64, @floatFromInt(lat));

            self.* = .{
                .src_sample_rate = src_rate,
                .dst_sample_rate = dst_rate,
                .init_frac_pos = init_frac,
                .latency = lat,
                .latency_frac = 0.0,
                .filter_len = 0,
                .fll = 0,
                .fl2 = 0,
                .flo = 0,
                .flb = 0,
                .filter_bank = undefined,
                .bank_cache = bank_cache,
                .is_whole = false,
                .is_bank_owned = false,
                .allocator = allocator,
            };

            var in_step: i32 = 0;
            var out_step: i32 = 0;
            if (getWholeStepping(src_rate, dst_rate, &in_step, &out_step)) {
                self.is_whole = true;
                self.in_step = in_step;
                self.out_step = out_step;
                const spos_w = self.init_frac_pos * @as(f64, @floatFromInt(out_step));
                self.init_frac_pos_w = @as(i32, @intFromFloat(spos_w));
                self.latency_frac = (spos_w - @as(f64, @floatFromInt(self.init_frac_pos_w))) / @as(f64, @floatFromInt(in_step));

                if (bank_cache) |bc| {
                    self.filter_bank = try bc.get(out_step, 1, 2, req_atten, is_third, correction_gain);
                } else {
                    self.filter_bank = try FracDelayFilterBank.init(allocator, out_step, 1, 2, req_atten, is_third, correction_gain);
                    self.is_bank_owned = true;
                }
            } else {
                self.is_whole = false;
                self.latency_frac = 0.0;
                self.frac_step = src_rate / dst_rate;

                if (bank_cache) |bc| {
                    self.filter_bank = try bc.get(-1, 3, 8, req_atten, is_third, correction_gain);
                } else {
                    self.filter_bank = try FracDelayFilterBank.init(allocator, -1, 3, 8, req_atten, is_third, correction_gain);
                    self.is_bank_owned = true;
                }
            }

            self.filter_len = self.filter_bank.filter_len;
            self.fl2 = self.filter_len >> 1;
            self.fll = self.fl2 - 1;
            self.flo = self.fll + self.fl2;
            self.flb = buf_len - self.fll;

            self.clear();
            return self;
        }

        pub fn deinit(self: *Self) void {
            if (self.bank_cache) |bc| {
                bc.unref(self.filter_bank);
            } else if (self.is_bank_owned) {
                self.filter_bank.deinit();
            }
            self.allocator.destroy(self);
        }

        pub fn getLatency(self: *Self) usize {
            return self.latency;
        }

        /// Unlike C++ (which omits fl2), this counts the fl2 samples convolve
        /// waits for before its first output — verified empirically.
        pub fn getInLenBeforeOutPos(self: *Self, req_out_pos: usize) usize {
            const ilat = self.fl2 + self.latency;
            if (self.is_whole) {
                return ilat + @as(usize, @intFromFloat((@as(f64, @floatFromInt(self.init_frac_pos_w)) + @as(f64, @floatFromInt(req_out_pos)) * @as(f64, @floatFromInt(self.in_step))) / @as(f64, @floatFromInt(self.out_step))));
            }
            return ilat + @as(usize, @intFromFloat(self.init_frac_pos + @as(f64, @floatFromInt(req_out_pos)) * self.src_sample_rate / self.dst_sample_rate));
        }

        pub fn getLatencyFrac(self: *Self) f64 {
            return self.latency_frac;
        }

        pub fn getName(ptr: *anyopaque) []const u8 {
            const self: *Self = @ptrCast(@alignCast(ptr));
            return if (self.is_whole) "FracInterpolator(Whole)" else "FracInterpolator(Fractional)";
        }

        pub fn getMaxOutLen(self: *Self, max_in_len: usize) usize {
            return @as(usize, @intFromFloat(@ceil(@as(f64, @floatFromInt(max_in_len)) * self.dst_sample_rate / self.src_sample_rate))) + 1;
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
            const out_len = self.getMaxOutLen(l);
            const written = self.process(ip[0..l], op0.*[0..out_len]);
            return written;
        }

        fn deinitGeneric(ptr: *anyopaque, allocator: std.mem.Allocator) void {
            _ = allocator;
            const self: *Self = @ptrCast(@alignCast(ptr));
            self.deinit();
        }

        pub fn clear(self: *Self) void {
            self.latency_left = self.latency;
            self.buf_left = 0;
            self.write_pos = 0;
            self.read_pos = self.flb;
            @memset(&self.buf, 0);

            if (self.is_whole) {
                self.in_pos_frac_w = self.init_frac_pos_w;
            } else {
                self.in_pos_frac = self.init_frac_pos;
                self.in_counter = 0;
                self.in_pos_int = 0;
                self.in_pos_shift = self.init_frac_pos * self.dst_sample_rate / self.src_sample_rate;
            }
        }

        fn convolve0(self: *Self, op_in: *[]T) void {
            // ponytail: comptime-dispatch on the one tap count per bank.
            switch (self.filter_len) {
                inline 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30 => |n| self.convolve0N(n, op_in),
                else => unreachable,
            }
        }

        fn convolve0N(self: *Self, comptime fltlen: usize, op_in: *[]T) void {
            const fb = self.filter_bank;
            const istep = self.in_step;
            const ostep = self.out_step;
            var fpos = self.in_pos_frac_w;
            var rpos = self.read_pos;
            var bl = @as(isize, @intCast(self.buf_left)) - @as(isize, @intCast(self.fl2));
            var op = op_in.*;

            while (bl > 0) {
                const ftp = fb.getFilter(@as(usize, @intCast(fpos)));
                op[0] = base.firDot(T, ftp[0..fltlen], self.buf[rpos..]);
                op = op[1..];

                fpos += istep;
                const pos_incr = @divTrunc(fpos, ostep);
                fpos -= pos_incr * ostep;

                rpos = (rpos + @as(usize, @intCast(pos_incr))) & buf_len_mask;
                bl -= pos_incr;
                if (op.len == 0) break;
            }

            self.buf_left = @intCast(bl + @as(isize, @intCast(self.fl2)));
            self.read_pos = rpos;
            self.in_pos_frac_w = fpos;
            op_in.* = op;
        }

        /// Σ (f0 + f1*x + f2*x²)[i] * rp[i]; ftp is planar per fir_lanes block (see FracDelayFilterBank.init).
        inline fn convolveQuadratic(ftp: []const f64, rp: []const T, comptime fltlen: usize, x: f64, x2: f64) T {
            const L = base.fir_lanes;
            const V = @Vector(L, f64);
            const xv: V = @splat(x);
            const x2v: V = @splat(x2);
            var acc: @Vector(L, T) = @splat(0);
            var i: usize = 0;
            while (i + L <= fltlen) : (i += L) {
                const f0: V = ftp[i * 3 ..][0..L].*;
                const f1: V = ftp[i * 3 + L ..][0..L].*;
                const f2: V = ftp[i * 3 + 2 * L ..][0..L].*;
                const f: @Vector(L, T) = @floatCast(f0 + f1 * xv + f2 * x2v);
                acc += f * @as(@Vector(L, T), rp[i..][0..L].*);
            }
            var s = @reduce(.Add, acc);
            while (i < fltlen) : (i += 1) {
                s += @as(T, @floatCast(ftp[i * 3] + ftp[i * 3 + 1] * x + ftp[i * 3 + 2] * x2)) * rp[i];
            }
            return s;
        }

        fn convolve2(self: *Self, op_in: *[]T) void {
            switch (self.filter_len) {
                inline 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30 => |n| self.convolve2N(n, op_in),
                else => unreachable,
            }
        }

        fn convolve2N(self: *Self, comptime fltlen: usize, op_in: *[]T) void {
            const fb = self.filter_bank;
            const ssr = self.src_sample_rate;
            const dsr = self.dst_sample_rate;
            var fpos = self.in_pos_frac;
            var rpos = self.read_pos;
            var bl = @as(isize, @intCast(self.buf_left)) - @as(isize, @intCast(self.fl2));
            var op = op_in.*;

            const fracs_f = @as(f64, @floatFromInt(fb.filter_fracs));

            while (bl > 0) {
                var x = fpos * fracs_f;
                const fti = @as(usize, @intFromFloat(x));
                x -= @as(f64, @floatFromInt(fti));
                const x2d = x * x;
                
                const ftp = fb.getFilter(fti);
                op[0] = convolveQuadratic(ftp, self.buf[rpos..], fltlen, x, x2d);
                op = op[1..];

                self.in_counter += 1;
                const next_in_pos = (@as(f64, @floatFromInt(self.in_counter)) + self.in_pos_shift) * ssr / dsr;
                const next_in_pos_int = @as(i32, @intFromFloat(@floor(next_in_pos)));
                const pos_incr = next_in_pos_int - self.in_pos_int;
                self.in_pos_int = next_in_pos_int;
                fpos = next_in_pos - @as(f64, @floatFromInt(next_in_pos_int));

                rpos = (rpos + @as(usize, @intCast(pos_incr))) & buf_len_mask;
                bl -= pos_incr;
                if (op.len == 0) break;
            }

            self.buf_left = @intCast(bl + @as(isize, @intCast(self.fl2)));
            self.read_pos = rpos;
            self.in_pos_frac = fpos;
            op_in.* = op;
        }

        pub fn process(self: *Self, ip_in: []const T, op_full: []T) usize {
            var l = ip_in.len;
            var ip = ip_in;
            var op = op_full;

            if (self.latency_left != 0) {
                const skip = @min(self.latency_left, l);
                self.latency_left -= skip;
                ip = ip[skip..];
                l -= skip;
            }

            while (l > 0) {
                const b = @min(l, @min(buf_len - self.write_pos, self.flb - self.buf_left));
                @memcpy(self.buf[self.write_pos .. self.write_pos + b], ip[0..b]);
                
                const ec = @as(isize, @intCast(self.flo)) - @as(isize, @intCast(self.write_pos));
                if (ec > 0) {
                    const copy_len = @min(b, @as(usize, @intCast(ec)));
                    @memcpy(self.buf[self.write_pos + buf_len .. self.write_pos + buf_len + copy_len], ip[0..copy_len]);
                }

                ip = ip[b..];
                self.write_pos = (self.write_pos + b) & buf_len_mask;
                l -= b;
                self.buf_left += b;

                if (self.is_whole) {
                    self.convolve0(&op);
                } else {
                    self.convolve2(&op);
                }
            }

            // C++: rebase the timing counter for precision (and so the i32
            // never overflows on long streams).
            if (!self.is_whole and self.in_counter > 1000) {
                self.in_counter = 0;
                self.in_pos_int = 0;
                self.in_pos_shift = self.in_pos_frac * self.dst_sample_rate / self.src_sample_rate;
            }

            const ol = op_full.len - op.len;
            return ol;
        }
    };
}

fn findGCD(l_in: f64, s_in: f64, gcd: *f64) bool {
    var l = l_in;
    var s = s_in;
    var it: i32 = 0;
    while (it < 150) : (it += 1) {
        const r = l - s;
        if (r == 0.0) {
            gcd.* = s;
            return s > 0.0;
        }
        l = s;
        s = @abs(r);
    }
    return false;
}

pub fn getWholeStepping(src_rate: f64, dst_rate: f64, in_step: *i32, out_step: *i32) bool {
    var gcd: f64 = 0;
    if (!findGCD(src_rate, dst_rate, &gcd)) return false;

    const in_step0 = src_rate / gcd;
    in_step.* = @as(i32, @intFromFloat(in_step0));
    const out_step0 = dst_rate / gcd;
    out_step.* = @as(i32, @intFromFloat(out_step0));

    if (@as(f64, @floatFromInt(in_step.*)) != in_step0 or @as(f64, @floatFromInt(out_step.*)) != out_step0) {
        return false;
    }

    if (out_step.* > 1500) return false;
    return true;
}
