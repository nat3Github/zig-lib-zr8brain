const std = @import("std");
const base = @import("base.zig");
const filters = @import("filters.zig");
const halfband = @import("halfband.zig");
const interpolator = @import("interpolator.zig");

pub fn Resampler(comptime T: type) type {
    return struct {
        const Self = @This();

        processors: std.ArrayList(base.Processor),
        tmp_bufs: [2][]f64,
        tmp_buf_all: []f64,
        tmp_result: []T,
        max_in_len: usize,
        latency_frac: f64,
        total_latency: usize,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator, src_rate: f64, dst_rate: f64, max_in_len: usize, steepness: f64, req_atten: f64, req_phase: filters.FilterPhaseResponse, cache: *filters.FIRFilterCache) !*Self {
            const self = try allocator.create(Self);
            self.* = .{
                .processors = .empty,
                .tmp_bufs = .{ &.{}, &.{} },
                .tmp_buf_all = &.{},
                .tmp_result = &.{},
                .max_in_len = max_in_len,
                .latency_frac = 0.0,
                .total_latency = 0,
                .allocator = allocator,
            };
            errdefer self.deinit();

            if (src_rate == dst_rate) return try self.finalizeInit();

            // Match C++ logic: prioritize power-of-2 / half-band stages before common ratios

            var cur_lat: f64 = 0.0;
            // Power-of-2 or 3*power-of-2 upsampling
            var i: usize = 2;
            while (i <= 3) : (i += 1) {
                var was_found = false;
                var c: usize = 0;
                while (true) {
                    const new_sr = src_rate * @as(f64, @floatFromInt(i << @as(std.math.Log2Int(usize), @intCast(c))));
                    if (new_sr == dst_rate) {
                        was_found = true;
                        break;
                    }
                    if (new_sr > dst_rate) break;
                    c += 1;
                }

                if (was_found) {
                    cur_lat = 0.0;
                    const filter = try cache.getLPFilter(1.0 / @as(f64, @floatFromInt(i)), steepness, req_atten, req_phase, @as(f64, @floatFromInt(i)));
                    const conv = try filters.BlockConvolver(f64).init(allocator, filter, @as(u32, @intCast(i)), 1, cur_lat, true, cache);
                    const proc = conv.processor();
                    try self.addProcessor(proc);
                    cur_lat = proc.getLatencyFrac();

                    const is_third = (i == 3);
                    var k: usize = 0;
                    while (k < c) : (k += 1) {
                        const up = try allocator.create(halfband.HBUpsampler(f64));
                        up.* = halfband.HBUpsampler(f64).init(req_atten, @as(i32, @intCast(k)), is_third, cur_lat, true);
                        const proc_up = up.processor();
                        try self.addProcessor(proc_up);
                        cur_lat = proc_up.getLatencyFrac();
                    }
                    return try self.finalizeInit();
                }
            }

            if (dst_rate * 2.0 > src_rate) {
                cur_lat = 0.0;
                // Upsampling or fractional downsampling down to 2X
                const norm_freq = if (dst_rate > src_rate) 0.5 else 0.5 * dst_rate / src_rate;
                const filter = try cache.getLPFilter(norm_freq, steepness, req_atten, req_phase, 2.0);
                const conv = try filters.BlockConvolver(f64).init(allocator, filter, 2, 1, cur_lat, true, cache);
                const proc = conv.processor();
                try self.addProcessor(proc);
                cur_lat = proc.getLatencyFrac();

                const tbw = 0.0175;
                const thresh_rate = src_rate / (1.0 - tbw * steepness);

                var c: usize = 0;
                var div: usize = 1;
                while (true) {
                    const ndiv = div * 2;
                    if (dst_rate < thresh_rate * @as(f64, @floatFromInt(ndiv))) break;
                    div = ndiv;
                    c += 1;
                }

                var c2: usize = 0;
                var div2: usize = 1;
                while (true) {
                    const ndiv = div * (if (c2 == 0) @as(usize, 3) else @as(usize, 2));
                    if (dst_rate < thresh_rate * @as(f64, @floatFromInt(ndiv))) break;
                    div2 = ndiv;
                    c2 += 1;
                }

                const src_rate2 = src_rate * 2.0;
                var tmp1: i32 = 0;
                var tmp2: i32 = 0;
                if (c == 1 and interpolator.getWholeStepping(src_rate2, dst_rate, &tmp1, &tmp2)) {
                    c = 0;
                }

                if (c > 0) {
                    // Steps using intermediate interpolation — faithful C++
                    // port of CDSPResampler: FracInterpolator → BlockConvolver
                    // (↑num) → (c-1)× HBUpsampler. The earlier port dropped the
                    // BlockConvolver stage entirely and ran the others in the
                    // wrong order, producing garbage for every c>0 ratio.
                    var num: usize = undefined;
                    if (c2 > 0 and div2 > div) {
                        div = div2;
                        c = c2;
                        num = 3;
                    } else {
                        num = 2;
                    }
                    const is_third = (num == 3);

                    const interp = try interpolator.FracInterpolator(f64).initWithCache(
                        allocator,
                        src_rate2 * @as(f64, @floatFromInt(div)),
                        dst_rate,
                        req_atten,
                        false,
                        cur_lat,
                        1.0,
                        &cache.bank_cache,
                    );
                    const interp_proc = interp.processor();
                    try self.addProcessor(interp_proc);
                    cur_lat = interp_proc.getLatencyFrac();

                    // Intermediate filter's transition band — kept no steeper
                    // than ReqTransBand, clamped to getLPMaxTransBand() (45.0).
                    var tb = (1.0 - src_rate * @as(f64, @floatFromInt(div)) / dst_rate) / tbw;
                    if (tb > 45.0) tb = 45.0;
                    const conv_filter = try cache.getLPFilter(
                        1.0 / @as(f64, @floatFromInt(num)),
                        tb,
                        req_atten,
                        req_phase,
                        @floatFromInt(num),
                    );
                    const interp_conv = try filters.BlockConvolver(f64).init(
                        allocator,
                        conv_filter,
                        @intCast(num),
                        1,
                        cur_lat,
                        true,
                        cache,
                    );
                    const interp_conv_proc = interp_conv.processor();
                    try self.addProcessor(interp_conv_proc);
                    cur_lat = interp_conv_proc.getLatencyFrac();

                    var hi: usize = 1;
                    while (hi < c) : (hi += 1) {
                        const up = try allocator.create(halfband.HBUpsampler(f64));
                        up.* = halfband.HBUpsampler(f64).init(req_atten, @as(i32, @intCast(hi - 1)), is_third, cur_lat, true);
                        const proc_up = up.processor();
                        try self.addProcessor(proc_up);
                        cur_lat = proc_up.getLatencyFrac();
                    }
                } else {
                    const interp = try interpolator.FracInterpolator(f64).initWithCache(allocator, src_rate2, dst_rate, req_atten, false, cur_lat, 1.0, &cache.bank_cache);
                    try self.addProcessor(interp.processor());
                }
                return try self.finalizeInit();
            } else {
                // Downsampling pipeline — faithful C++ port of CDSPResampler downsampling branch.
                // C++ threshold: HBDownsampler stages are added only when src >= 4 * dst * 2^c.
                var check_sr = dst_rate * 4.0;
                var c: usize = 0;
                var fin_gain: f64 = 1.0;

                while (check_sr <= src_rate) {
                    c += 1;
                    check_sr *= 2.0;
                    fin_gain *= 0.5;
                }

                const src_sr_div: usize = @as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(c));

                // Determine integer downsampling factor and whether fractional interpolation is needed.
                var downf: usize = 1;
                var norm_freq: f64 = 0.5;
                var use_interp = true;
                var is_third = false;

                var df: usize = 2;
                while (df <= 3) : (df += 1) {
                    if (dst_rate * @as(f64, @floatFromInt(src_sr_div)) * @as(f64, @floatFromInt(df)) == src_rate) {
                        norm_freq = 1.0 / @as(f64, @floatFromInt(df));
                        use_interp = false;
                        is_third = (df == 3);
                        downf = df;
                        break;
                    }
                }

                if (use_interp) {
                    norm_freq = dst_rate * @as(f64, @floatFromInt(src_sr_div)) / src_rate;
                    is_third = (norm_freq * 3.0 <= 1.0);
                }

                // Add HBDownsampler stages (each halves the sample rate).
                cur_lat = 0.0;
                var si: usize = 0;
                while (si < c) : (si += 1) {
                    const k = c - 1 - si;
                    const down = try allocator.create(halfband.HBDownsampler(f64));
                    down.* = halfband.HBDownsampler(f64).init(req_atten, @as(i32, @intCast(k)), is_third, cur_lat);
                    const proc = down.processor();
                    try self.addProcessor(proc);
                    cur_lat = proc.getLatencyFrac();
                }

                // Add LP BlockConvolver (with optional integer downsampling by downf).
                const filter = try cache.getLPFilter(norm_freq, steepness, req_atten, req_phase, fin_gain);
                const conv = try filters.BlockConvolver(f64).init(allocator, filter, 1, @as(u32, @intCast(downf)), cur_lat, true, cache);
                const conv_proc = conv.processor();
                try self.addProcessor(conv_proc);
                cur_lat = conv_proc.getLatencyFrac();

                // Add FracInterpolator for any remaining fractional ratio.
                if (use_interp) {
                    const interp = try interpolator.FracInterpolator(f64).initWithCache(
                        allocator,
                        src_rate,
                        dst_rate * @as(f64, @floatFromInt(src_sr_div)),
                        req_atten, is_third, cur_lat, 1.0, &cache.bank_cache);
                    try self.addProcessor(interp.processor());
                }
                return try self.finalizeInit();
            }
        }

        fn finalizeInit(self: *Self) !*Self {
            try self.createTmpBuffers();
            return self;
        }

        pub fn deinit(self: *Self) void {
            for (self.processors.items) |p| {
                p.deinit(self.allocator);
            }
            self.processors.deinit(self.allocator);
            if (self.tmp_buf_all.len > 0) self.allocator.free(self.tmp_buf_all);
            if (self.tmp_result.len > 0) self.allocator.free(self.tmp_result);
            self.allocator.destroy(self);
        }

        /// Takes ownership of `proc`, also on failure.
        fn addProcessor(self: *Self, proc: base.Processor) !void {
            self.processors.append(self.allocator, proc) catch |err| {
                proc.deinit(self.allocator);
                return err;
            };
            self.latency_frac = proc.getLatencyFrac();
        }

        fn createTmpBuffers(self: *Self) !void {
            var cur_max_out = self.max_in_len;
            // tmp_bufs[0] holds the initial input and output of odd-indexed stages.
            // tmp_bufs[1] holds the output of even-indexed stages (0, 2, ...).
            // Stage i writes to tmp_bufs[(i+1)&1], so stage 0 → tmp_bufs[1].
            var caps = [2]usize{ self.max_in_len, 0 };
            var cur_buf: usize = 1; // stage 0's output goes to tmp_bufs[1]

            for (self.processors.items) |p| {
                cur_max_out = p.getMaxOutLen(cur_max_out);
                if (cur_max_out > caps[cur_buf]) {
                    caps[cur_buf] = cur_max_out;
                }
                cur_buf ^= 1;
            }

            const total = caps[0] + caps[1];
            if (total > 0) {
                self.tmp_buf_all = try self.allocator.alloc(f64, total);
                @memset(self.tmp_buf_all, 0);
                self.tmp_bufs[0] = self.tmp_buf_all[0..caps[0]];
                self.tmp_bufs[1] = self.tmp_buf_all[caps[0]..];
                
                // tmp_result for type conversion if needed
                self.tmp_result = try self.allocator.alloc(T, cur_max_out);
                @memset(self.tmp_result, 0);
            }
        }

        pub fn process(self: *Self, ip0: []const T, op_full: []T) usize {
            // Boundary NaN check: caller bug if NaN reaches the pipeline.
            // Stripped in ReleaseFast/ReleaseSmall.
            if (std.debug.runtime_safety) {
                for (ip0) |s| std.debug.assert(!std.math.isNan(s));
            }

            const out = if (T == f64) self.processRaw(ip0) else blk: {
                const in64 = self.tmp_bufs[0][0..ip0.len];
                for (in64, ip0) |*d, s| d.* = @floatCast(s);
                break :blk self.processRaw(in64);
            };
            const out_len = @min(out.len, op_full.len);
            for (op_full[0..out_len], out[0..out_len]) |*d, s| d.* = @floatCast(s);
            return out_len;
        }

        /// C++ CDSPResampler::process semantics: `ip` is read in place and the
        /// result aliases an internal buffer, valid until the next call.
        pub fn processRaw(self: *Self, ip0: []const f64) []const f64 {
            var l = ip0.len;
            var ip = ip0.ptr;
            for (self.processors.items, 0..) |p, i| {
                const op_start = self.tmp_bufs[(i + 1) & 1].ptr;
                var op_ptr = op_start;
                l = p.process(ip, l, &op_ptr);
                ip = op_start;
            }
            return ip[0..l];
        }

        pub fn clear(self: *Self) void {
            for (self.processors.items) |p| {
                p.clear();
            }
        }

        pub fn debugPrintPipeline(self: *Self) void {
            std.debug.print("--- Resampler Pipeline Trace ---\n", .{});
            for (self.processors.items, 0..) |p, i| {
                const name = p.getName();
                if (std.mem.eql(u8, name, "BlockConvolver")) {
                    const bc: *filters.BlockConvolver(T) = @ptrCast(@alignCast(p.ptr));
                    std.debug.print("  Stage[{d}]: {s} | Ptr: {*} | LatencyFrac: {d:.10} | InLen(0): {d} | Block: {d} | InLen: {d} | Delay: {d}\n", 
                        .{ i, name, p.ptr, p.getLatencyFrac(), p.getInLenBeforeOutPos(0), bc.block_len2, bc.input_len, bc.input_delay });
                } else {
                    std.debug.print("  Stage[{d}]: {s} | Ptr: {*} | LatencyFrac: {d:.10} | InLen(0): {d}\n", .{ i, name, p.ptr, p.getLatencyFrac(), p.getInLenBeforeOutPos(0) });
                }
            }
            std.debug.print("  Total Reported Latency (Input scale, pos 0): {d}\n", .{self.getInLenBeforeOutPos(0)});
            std.debug.print("--------------------------------\n", .{});
        }

        pub fn getLatency(self: *Self) usize {
            return self.getInLenBeforeOutPos(0);
        }

        pub fn getInLenBeforeOutPos(self: *Self, req_out_pos: usize) usize {
            var req = req_out_pos;
            var i: usize = self.processors.items.len;
            while (i > 0) {
                i -= 1;
                req = self.processors.items[i].getInLenBeforeOutPos(req);
            }
            return req;
        }

        pub fn getInputRequiredForOutput(self: *Self, req_out_samples: usize) usize {
            if (req_out_samples < 1) return 0;
            return self.getInLenBeforeOutPos(req_out_samples - 1) + 1;
        }

        pub fn getMaxOutLen(self: *Self, max_in: usize) usize {
            var l = max_in;
            for (self.processors.items) |p| {
                l = p.getMaxOutLen(l);
            }
            return l;
        }
    };
}
