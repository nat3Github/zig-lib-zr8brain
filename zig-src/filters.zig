const std = @import("std");
const base = @import("base.zig");
const fft = @import("fft.zig");

pub const FilterPhaseResponse = enum {
    linearPhase,
    minPhase,
};

pub fn SincFilterGen(comptime T: type) type {
    return struct {
        const Self = @This();
        
        len2: T,
        len2i: T = 0,
        kernel_len: usize = 0,
        fl2: usize = 0,
        
        freq1: T = 0,
        freq2: T = 0,
        frac_delay: T = 0,
        
        power: T = 0,
        w1: base.SineGen = undefined,
        w2: base.SineGen = undefined,
        w3: base.SineGen = undefined,
        
        window_type: WindowFunctionType = .cosine,
        
        kaiser_beta: T = 0,
        kaiser_mul: T = 0,
        kaiser_len2_frac: T = 0,
        
        gaussian_sigma_i: T = 0,
        gaussian_sigma_frac: T = 0,
        
        wn: isize = 0,

        pub const WindowFunctionType = enum {
            cosine,
            kaiser,
            gaussian,
        };

        pub fn calcWindowKaiser(self: *Self) T {
            const x = self.kaiser_len2_frac + @as(T, @floatFromInt(self.wn)) * self.len2i;
            const arg = 1.0 - x * x;
            const win = besselI0(self.kaiser_beta * std.math.sqrt(@max(0.0, arg))) * self.kaiser_mul;
            // arg >= 0 via @max above; besselI0 finite for finite input → win finite
            if (std.math.isNan(win)) unreachable;
            self.wn += 1;
            return win;
        }

        fn besselI0(x: T) T {
            const ax = @abs(x);
            if (ax < 3.75) {
                const y = x / 3.75;
                const ys = y * y;
                return 1.0 + ys * (3.5156229 + ys * (3.0899424 + ys * (1.2067492 + ys * (0.2659732 + ys * (0.0360768 + ys * 0.0045813)))));
            }
            const y = 3.75 / ax;
            return @exp(ax) / @sqrt(ax) * (0.39894228 + y * (0.01328592 + y * (0.00225319 + y * (-0.00157565 + y * (0.00916281 + y * (-0.02057706 + y * (0.02635537 + y * (-0.01647633 + y * 0.00392377))))))));
        }

        pub fn initBand(self: *Self, win_type: WindowFunctionType, params: ?[]const T, use_power: bool) void {
            std.debug.assert(self.len2 >= 2.0);
            self.fl2 = @intFromFloat(@floor(self.len2));
            self.kernel_len = self.fl2 + self.fl2 + 1;
            self.setWindow(win_type, params, use_power, true, 0.0);
        }

        pub fn initFrac(self: *Self, win_type: WindowFunctionType, params: ?[]const T, use_power: bool, kernel_len: usize) void {
            std.debug.assert(self.len2 >= 2.0);
            self.fl2 = @intFromFloat(@floor(self.len2));
            self.kernel_len = kernel_len;
            self.setWindow(win_type, params, use_power, false, self.frac_delay);
        }

        fn setWindow(self: *Self, win_type: WindowFunctionType, params: ?[]const T, use_power: bool, is_centered: bool, use_frac_delay: T) void {
            self.frac_delay = use_frac_delay;
            self.window_type = win_type;

            if (win_type == .cosine) {
                @panic("TODO");
            } else if (win_type == .kaiser) {
                self.wn = if (is_centered) 0 else -@as(isize, @intCast(self.fl2));
                if (params == null) {
                    self.kaiser_beta = 9.5945013206755156;
                    self.power = if (use_power) 1.9718457932433306 else -1.0;
                } else {
                    self.kaiser_beta = std.math.clamp(params.?[0], 1.0, 350.0);
                    self.power = if (use_power) @abs(params.?[1]) else -1.0;
                }
                self.kaiser_mul = 1.0 / besselI0(self.kaiser_beta);
                self.len2i = 1.0 / self.len2;
                self.kaiser_len2_frac = self.frac_delay * self.len2i;
            } else if (win_type == .gaussian) {
                @panic("TODO");
            }
        }

        pub fn generateBand(self: *Self, op_full: []T) void {
            const pi = std.math.pi;
            var f2 = base.SineGen.init(self.freq2, 0.0, 1.0 / pi);
            _ = f2.generate();

            const c_idx = self.fl2;
            const op = op_full;
            const pw = self.power;
            var t: usize = 1;

            if (self.freq1 < 2.3e-13) {
                if (pw < 0.0) {
                    op[c_idx] = self.freq2 * self.calcWindowKaiser() / pi;
                    while (t <= self.fl2) : (t += 1) {
                        const v = f2.generate() * self.calcWindowKaiser() / @as(T, @floatFromInt(t));
                        op[c_idx + t] = v;
                        op[c_idx - t] = v;
                    }
                } else {
                    op[c_idx] = self.freq2 * std.math.pow(T, self.calcWindowKaiser(), pw) / pi;
                    while (t <= self.fl2) : (t += 1) {
                        const v = f2.generate() * std.math.pow(T, self.calcWindowKaiser(), pw) / @as(T, @floatFromInt(t));
                        op[c_idx + t] = v;
                        op[c_idx - t] = v;
                    }
                }
            } else {
                var f1 = base.SineGen.init(self.freq1, 0.0, 1.0 / pi);
                _ = f1.generate();

                if (pw < 0.0) {
                    op[c_idx] = (self.freq2 - self.freq1) * self.calcWindowKaiser() / pi;
                    while (t <= self.fl2) : (t += 1) {
                        const v = (f2.generate() - f1.generate()) * self.calcWindowKaiser() / @as(T, @floatFromInt(t));
                        op[c_idx + t] = v;
                        op[c_idx - t] = v;
                    }
                } else {
                    op[c_idx] = (self.freq2 - self.freq1) * std.math.pow(T, self.calcWindowKaiser(), pw) / pi;
                    while (t <= self.fl2) : (t += 1) {
                        const v = (f2.generate() - f1.generate()) * std.math.pow(T, self.calcWindowKaiser(), pw) / @as(T, @floatFromInt(t));
                        op[c_idx + t] = v;
                        op[c_idx - t] = v;
                    }
                }
            }
        }

        pub fn generateFrac(self: *Self, op: []T, element_size: usize) void {
            const pi = std.math.pi;
            const f2 = self.freq2;
            const f1 = self.freq1;
            const pw = self.power;
            var t_val = -self.len2 + self.frac_delay;

            for (0..self.kernel_len) |i| {
                var v: T = undefined;
                if (t_val == 0.0) {
                    v = f2 / pi;
                } else {
                    v = @sin(f2 * t_val) / (pi * t_val);
                }

                if (f1 > 0.0) {
                    if (t_val == 0.0) {
                        v -= f1 / pi;
                    } else {
                        v -= @sin(f1 * t_val) / (pi * t_val);
                    }
                }

                const win = if (pw < 0.0) self.calcWindowKaiser() else std.math.pow(T, self.calcWindowKaiser(), pw);
                op[i * element_size] = v * win;
                t_val += 1.0;
            }
        }
    };
}

fn calcFIRFilterResponse(flt: []const f64, th: f64) [2]f64 {
    const sincr = 2.0 * @cos(th);
    var c1: f64 = 1.0;
    var s1: f64 = 0.0;
    var c2: f64 = @cos(-th);
    var s2: f64 = @sin(-th);
    var re: f64 = 0.0;
    var im: f64 = 0.0;
    for (flt) |v| {
        re += c1 * v;
        im += s1 * v;
        const tc = c1;
        c1 = sincr * c1 - c2;
        c2 = tc;
        const ts = s1;
        s1 = sincr * s1 - s2;
        s2 = ts;
    }
    return .{ re, im };
}

fn calcFIRFilterGroupDelay(flt: []const f64, th: f64) f64 {
    const thd2 = 1e-9;
    var ths0 = th - thd2;
    var ths1 = th + thd2;
    if (ths0 < 0.0) ths0 = 0.0;
    if (ths1 > std.math.pi) ths1 = std.math.pi;
    const r0 = calcFIRFilterResponse(flt, ths0);
    const r1 = calcFIRFilterResponse(flt, ths1);
    const ph0 = std.math.atan2(r0[1], r0[0]);
    var ph1 = std.math.atan2(r1[1], r1[0]);
    if (@abs(ph1 - ph0) > std.math.pi) {
        if (ph1 > ph0) ph1 -= 2.0 * std.math.pi else ph1 += 2.0 * std.math.pi;
    }
    return (ph1 - ph0) / (ths1 - ths0);
}

/// Cepstrum-based min-phase transform (mirrors CDSPRealFFT::calcMinPhaseTransform).
/// LenMult=16, DoFinalMul=false matches CDSPFIRFilter usage.
/// Returns DC group delay in samples.
fn calcMinPhaseTransform(allocator: std.mem.Allocator, kernel: []f64, len_mult: usize, do_final_mul: bool) !f64 {
    const kernel_len = kernel.len;
    const len_bits = base.getBitOccupancy(kernel_len * len_mult - 1);
    const len = @as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(len_bits));
    const len2 = len >> 1;

    var ip = try allocator.alloc(f64, len);
    defer allocator.free(ip);
    var ip2 = try allocator.alloc(f64, len2 + 1);
    defer allocator.free(ip2);

    @memcpy(ip[0..kernel_len], kernel);
    @memset(ip[kernel_len..len], 0);

    var ffto = try fft.RealFFT(f64).init(allocator, len_bits);
    defer ffto.deinit();
    const work = try ffto.allocWork(allocator);
    defer allocator.free(work);

    ffto.forward(ip, work);

    const nzbias = 1e-300;

    ip2[0] = ip[0];
    ip[0] = @log(@abs(ip[0]) + nzbias);
    ip2[len2] = ip[1];
    ip[1] = @log(@abs(ip[1]) + nzbias);

    for (1..len2) |i| {
        ip2[i] = @sqrt(ip[i * 2] * ip[i * 2] + ip[i * 2 + 1] * ip[i * 2 + 1]);
        ip[i * 2] = @log(ip2[i] + nzbias);
        ip[i * 2 + 1] = 0.0;
    }

    ffto.inverse(ip, work);

    const m1 = ffto.inv_mul_const;
    const m2 = -m1;

    ip[0] = 0.0;
    for (1..len2) |i| ip[i] *= m1;
    ip[len2] = 0.0;
    for ((len2 + 1)..len) |i| ip[i] *= m2;

    ffto.forward(ip, work);

    ip[0] = ip2[0];
    ip[1] = ip2[len2];
    for (1..len2) |i| {
        const phase = ip[i * 2 + 1];
        const mag = ip2[i];
        ip[i * 2] = @cos(phase) * mag;
        ip[i * 2 + 1] = @sin(phase) * mag;
    }

    ffto.inverse(ip, work);

    if (do_final_mul) {
        for (0..kernel_len) |i| kernel[i] = ip[i] * m1;
    } else {
        @memcpy(kernel, ip[0..kernel_len]);
    }

    return calcFIRFilterGroupDelay(kernel, 0.0);
}

pub fn FIRFilter(comptime T: type) type {
    return struct {
        const Self = @This();
        const RFFT = fft.RealFFT(T);

        req_norm_freq: T,
        req_trans_band: T,
        req_atten: T,
        req_phase: FilterPhaseResponse,
        req_gain: T,

        is_zero_phase: bool = false,
        latency: usize = 0,
        latency_frac: T = 0.0,
        kernel_len: usize = 0,
        block_len_bits: usize = 0,
        
        kernel_block: []align(64) T,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator, req_norm_freq: T, req_trans_band: T, req_atten: T, req_phase: FilterPhaseResponse, req_gain: T) !*Self {
            var self = try allocator.create(Self);
            self.* = .{
                .req_norm_freq = req_norm_freq,
                .req_trans_band = req_trans_band,
                .req_atten = req_atten,
                .req_phase = req_phase,
                .req_gain = req_gain,
                .kernel_block = &.{},
                .allocator = allocator,
            };

            errdefer self.deinit();
            try self.buildLPFilter(null);
            return self;
        }

        pub fn deinit(self: *Self) void {
            if (self.kernel_block.len > 0) {
                self.allocator.free(self.kernel_block);
            }
            self.allocator.destroy(self);
        }

        fn buildLPFilter(self: *Self, ext_atten_corrs: ?[]const T) !void {
            const tb = self.req_trans_band * 0.01;
            var atten = -self.req_atten;

            if (tb >= 0.25) {
                if (self.req_atten >= 117.0) { atten -= 1.60; } else if (self.req_atten >= 60.0) { atten -= 1.91; } else { atten -= 2.25; }
            } else if (tb >= 0.10) {
                if (self.req_atten >= 117.0) { atten -= 0.69; } else if (self.req_atten >= 60.0) { atten -= 0.73; } else { atten -= 1.13; }
            } else {
                if (self.req_atten >= 117.0) { atten -= 0.21; } else if (self.req_atten >= 60.0) { atten -= 0.25; } else { atten -= 0.36; }
            }

            const atten_corr_count = 264;
            const atten_corr_min = 49.0;
            const atten_corr_diff = 176.25;
            var atten_corr_idx = @as(isize, @intFromFloat(@floor((-atten - atten_corr_min) * atten_corr_count / atten_corr_diff + 0.5)));
            atten_corr_idx = @max(0, @min(atten_corr_idx, atten_corr_count));

            if (ext_atten_corrs) |ext| {
                atten -= ext[@intCast(atten_corr_idx)];
            } else if (tb >= 0.25) {
                const atten_corrs = [_]i8{-127, -127, -125, -125, -122, -119, -115, -110, -104, -97, -91, -82, -75, -24, -16, -6, 4, 14, 24, 29, 30, 32, 37, 44, 51, 57, 63, 67, 65, 50, 53, 56, 58, 60, 63, 64, 66, 68, 74, 77, 78, 78, 78, 79, 79, 60, 60, 60, 61, 59, 52, 47, 41, 36, 30, 24, 17, 9, 0, -8, -10, -11, -14, -13, -18, -25, -31, -38, -44, -50, -57, -63, -68, -74, -81, -89, -96, -101, -104, -107, -109, -110, -86, -84, -85, -82, -80, -77, -73, -67, -62, -55, -48, -42, -35, -30, -20, -11, -2, 5, 6, 6, 7, 11, 16, 21, 26, 34, 41, 46, 49, 52, 55, 56, 48, 49, 51, 51, 52, 52, 52, 52, 52, 51, 51, 50, 47, 47, 50, 48, 46, 42, 38, 35, 31, 27, 24, 20, 16, 12, 11, 12, 10, 8, 4, -1, -6, -11, -16, -19, -17, -21, -24, -27, -32, -34, -37, -38, -40, -41, -40, -40, -42, -41, -44, -45, -43, -41, -34, -31, -28, -24, -21, -18, -14, -10, -5, -1, 2, 5, 8, 7, 4, 3, 2, 2, 4, 6, 8, 9, 9, 10, 10, 10, 10, 9, 8, 9, 11, 14, 13, 12, 11, 10, 8, 7, 6, 5, 3, 2, 2, -1, -1, -3, -3, -4, -4, -5, -4, -6, -7, -9, -5, -1, -1, 0, 1, 0, -2, -3, -4, -5, -5, -8, -13, -13, -13, -12, -13, -12, -11, -11, -9, -8, -7, -5, -3, -1, 2, 4, 6, 9, 10, 11, 14, 18, 21, 24, 27, 30, 34, 37, 37, 39, 40};
                atten -= @as(T, @floatFromInt(atten_corrs[@intCast(atten_corr_idx)])) / 101.0;
            } else if (tb >= 0.10) {
                const atten_corrs = [_]i8{-113, -118, -122, -125, -126, -97, -95, -92, -92, -89, -82, -75, -69, -48, -42, -36, -30, -22, -14, -5, -2, 1, 6, 13, 22, 28, 35, 41, 48, 55, 56, 56, 61, 65, 71, 77, 81, 83, 85, 85, 74, 74, 73, 72, 71, 70, 68, 64, 59, 56, 49, 52, 46, 42, 36, 32, 26, 20, 13, 7, -2, -6, -10, -15, -20, -27, -33, -38, -44, -43, -48, -53, -57, -63, -69, -73, -75, -79, -81, -74, -76, -77, -77, -78, -81, -80, -80, -78, -76, -65, -62, -59, -56, -51, -48, -44, -38, -33, -25, -19, -13, -5, -1, 2, 7, 13, 17, 21, 25, 30, 35, 40, 45, 50, 53, 56, 57, 55, 58, 59, 62, 64, 67, 67, 68, 68, 62, 61, 61, 59, 59, 57, 57, 55, 52, 48, 42, 38, 35, 31, 26, 20, 15, 13, 10, 7, 3, -2, -8, -13, -17, -23, -28, -34, -37, -40, -41, -45, -48, -50, -53, -57, -59, -62, -63, -63, -57, -57, -56, -56, -54, -54, -53, -49, -48, -41, -38, -33, -31, -26, -23, -18, -12, -9, -7, -7, -3, 0, 5, 9, 14, 16, 20, 22, 21, 23, 25, 27, 28, 29, 34, 33, 35, 33, 31, 30, 29, 29, 26, 26, 25, 24, 20, 19, 15, 10, 8, 4, 1, -2, -6, -10, -16, -19, -23, -26, -27, -30, -34, -39, -43, -47, -51, -52, -54, -56, -58, -59, -62, -63, -66, -65, -65, -64, -59, -57, -54, -52, -48, -44, -42, -37, -32, -22, -17, -10, -3, 5, 13, 22, 30, 40, 50, 60, 72};
                atten -= @as(T, @floatFromInt(atten_corrs[@intCast(atten_corr_idx)])) / 210.0;
            } else {
                const atten_corrs = [_]i8{-15, -17, -20, -20, -20, -21, -20, -16, -17, -18, -17, -13, -12, -11, -9, -7, -5, -4, -1, 1, 3, 4, 5, 6, 7, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 10, 11, 10, 10, 8, 10, 11, 10, 11, 11, 13, 14, 15, 19, 27, 26, 23, 18, 14, 8, 4, -2, -6, -12, -17, -23, -28, -33, -37, -42, -46, -49, -53, -57, -60, -61, -64, -65, -67, -66, -66, -66, -65, -64, -61, -59, -56, -52, -48, -42, -38, -31, -27, -19, -13, -7, -1, 8, 14, 22, 29, 37, 45, 52, 59, 66, 73, 80, 86, 91, 96, 100, 104, 108, 111, 114, 115, 117, 118, 120, 120, 118, 117, 114, 113, 111, 107, 103, 99, 95, 89, 84, 78, 72, 66, 60, 52, 44, 37, 30, 21, 14, 6, -3, -11, -18, -26, -34, -43, -51, -58, -65, -73, -78, -85, -90, -97, -102, -107, -113, -115, -118, -121, -125, -125, -126, -126, -126, -125, -124, -121, -119, -115, -111, -109, -101, -102, -95, -88, -81, -73, -67, -63, -54, -47, -40, -33, -26, -18, -11, -5, 2, 8, 14, 19, 25, 31, 36, 37, 43, 47, 49, 51, 52, 57, 57, 56, 57, 58, 58, 58, 57, 56, 52, 52, 50, 48, 44, 41, 39, 37, 33, 31, 26, 24, 21, 18, 14, 11, 8, 4, 2, -2, -5, -7, -9, -11, -13, -15, -16, -18, -19, -20, -23, -24, -24, -25, -27, -26, -27, -29, -30, -31, -32, -35, -36, -39, -40, -44, -46, -51, -54, -59, -63, -69, -76, -83, -91, -98};
                atten -= @as(T, @floatFromInt(atten_corrs[@intCast(atten_corr_idx)])) / 196.0;
            }

            const pwr = 7.43932822146293e-8 * (atten * atten) + 0.000102747434588003 * @cos(0.00785021930010397 * atten) * @cos(0.633854318781239 + 0.103208573657699 * atten) - 0.00798132247867036 - 0.000903555213543865 * atten - 0.0969365532127236 * @exp(0.0779275237937911 * atten) - 1.37304948662012e-5 * atten * @cos(0.00785021930010397 * atten);

            var hl: T = 0;
            var fo1: T = 0;

            const asinh = struct {
                fn eval(v: T) T { return @log(v + @sqrt(v * v + 1.0)); }
            }.eval;

            if (pwr <= 0.067665322581) {
                if (tb >= 0.25) {
                    hl = 2.6778150875894 / tb + 300.547590563091 * std.math.atan(std.math.atan(2.68959772209918 * pwr)) / (5.5099277187035 * tb - tb * std.math.tanh(@cos(asinh(atten))));
                    fo1 = 0.987205355829873 * tb + 1.00011788929851 * std.math.atan2(
                        -0.321432067051302 - 6.19131357321578 * @sqrt(pwr),
                        hl + -1.14861472207245 / (hl - 14.1821147585957) + std.math.pow(T, 0.9521145021664, std.math.pow(T, std.math.atan2(1.12018764830637, tb), 2.10988901686912 * hl - 20.9691278378345)));
                } else if (tb >= 0.10) {
                    hl = (1.56688617018066 + 142.064321294568 * pwr + 0.00419441117131136 * @cos(243.633511747297 * pwr) - 0.022953443903576 * atten - 0.026629568860284 * @cos(127.715550622571 * pwr)) / tb;
                    fo1 = 0.982299356642411 * tb + 0.999441744774215 * asinh((-0.361783054039583 - 5.80540593623676 * @sqrt(pwr)) / hl);
                } else {
                    hl = (2.45739657014937 + 269.183679500541 * pwr * @cos(5.73225668178813 + std.math.atan2(std.math.cosh(0.988861169868941 - 17.2201556280744 * pwr), 1.08340138240431 * pwr))) / tb;
                    fo1 = 2.291956939 * tb + 0.01942450693 * tb * tb * hl - 4.67538973161837 * pwr * tb - 1.668433124 * tb * std.math.pow(T, pwr, pwr);
                }
            } else {
                if (tb >= 0.25) {
                    hl = (1.50258368698213 + 158.556968859477 * asinh(pwr) * std.math.tanh(57.9466246871383 * std.math.tanh(pwr)) - 0.0105440479814834 * atten) / tb;
                    fo1 = 0.994024401639321 * tb + (-0.236282717577215 - 6.8724924545387 * @sqrt(@sin(pwr))) / hl;
                } else if (tb >= 0.10) {
                    hl = (1.50277377248945 + 158.222625721046 * asinh(pwr) * std.math.tanh(1.02875299001715 + 42.072277322604 * pwr) - 0.0108380943845632 * atten) / tb;
                    fo1 = 0.992539376734551 * tb + (-0.251747813037178 - 6.74159892452584 * @sqrt(std.math.tanh(std.math.tanh(@tan(pwr))))) / hl;
                } else {
                    hl = (1.15990238966306 * pwr - 5.02124037125213 * (pwr * pwr) - 0.158676856669827 * atten * @cos(1.1609073390614 * pwr - 6.33932586197475 * pwr * (pwr * pwr))) / tb;
                    fo1 = 0.867344453126885 * tb + 0.052693817907757 * tb * @log(pwr) + 0.0895511178735932 * tb * std.math.atan(59.7538527741309 * pwr) - 0.0745653568081453 * pwr * tb;
                }
            }
            hl = if (pwr <= 0.067665322581) (if (tb >= 0.25) (2.45739657014937 + 269.183679500541 * pwr * @cos(5.73225668178813 + std.math.atan2(std.math.cosh(1.111861169868941 - 17.2201556280744 * pwr), 1.08340138240431 * pwr))) / tb else (2.45739657014937 + 269.183679500541 * pwr * @cos(5.73225668178813 + std.math.atan2(std.math.cosh(0.988861169868941 - 17.2201556280744 * pwr), 1.08340138240431 * pwr))) / tb) else (if (tb >= 0.25) (1.50258368698213 + 158.556968859477 * std.math.asinh(pwr) * std.math.tanh(57.9466246871383 * std.math.tanh(pwr)) - 0.0105440479814834 * atten) / tb else if (tb >= 0.10) (1.50277377248945 + 158.222625721046 * std.math.asinh(pwr) * std.math.tanh(1.02875299001715 + 42.072277322604 * pwr) - 0.0108380943845632 * atten) / tb else (1.15990238966306 * pwr - 5.02124037125213 * (pwr * pwr) - 0.158676856669827 * atten * @cos(1.1609073390614 * pwr - 6.33932586197475 * pwr * (pwr * pwr))) / tb);

            var win_params: [2]T = .{ 125.0, pwr };

            var sinc = SincFilterGen(T){
                .len2 = 0.25 * hl / self.req_norm_freq,
                .freq1 = 0.0,
                .freq2 = std.math.pi * (1.0 - fo1) * self.req_norm_freq,
            };
            sinc.initBand(.kaiser, &win_params, true);

            self.kernel_len = sinc.kernel_len;
            
            // Extension amount per R8B_EXTFFT (we will assume 0 or 1 for simplicity of port). Usually 0 or 1.
            self.block_len_bits = base.getBitOccupancy(self.kernel_len - 1);
            const block_len = @as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(self.block_len_bits));

            self.kernel_block = try self.allocator.alignedAlloc(T, .@"64", block_len * 2);
            @memset(self.kernel_block, 0);
            
            sinc.generateBand(self.kernel_block[0..sinc.kernel_len]);

            if (self.req_phase == .linearPhase) {
                self.is_zero_phase = true;
                self.latency = sinc.fl2;
                self.latency_frac = 0.0;
            } else {
                self.is_zero_phase = false;
                const dc_group_delay = try calcMinPhaseTransform(self.allocator, self.kernel_block[0..self.kernel_len], 16, false);
                self.latency = @intFromFloat(dc_group_delay);
                self.latency_frac = @floatCast(dc_group_delay - @as(f64, @floatFromInt(self.latency)));
            }

            var ffto = try RFFT.init(self.allocator, self.block_len_bits + 1);
            defer ffto.deinit();
            const work = try ffto.allocWork(self.allocator);
            defer self.allocator.free(work);

            if (self.is_zero_phase) {
                var s: T = 0;
                for (self.kernel_block[0..self.kernel_len]) |v| {
                    s += v;
                }
                
                s = ffto.inv_mul_const * self.req_gain / s;

                // Time-shift the filter so that zero-phase response is produced
                for (0..sinc.fl2 + 1) |i| {
                    self.kernel_block[i] = self.kernel_block[sinc.fl2 + i] * s;
                }
                var i: usize = 1;
                while (i <= sinc.fl2) : (i += 1) {
                    self.kernel_block[block_len * 2 - i] = self.kernel_block[i];
                }
                @memset(self.kernel_block[sinc.fl2 + 1 .. block_len * 2 - sinc.fl2 ], 0);

                ffto.forward(self.kernel_block, work);
                ffto.convertToZP(self.kernel_block);
            } else {
                // Minimum phase: normalize kernel, zero-pad, forward FFT (no convertToZP)
                var dc_gain: T = 0;
                for (self.kernel_block[0..self.kernel_len]) |v| dc_gain += v;
                const scale = ffto.inv_mul_const * self.req_gain / dc_gain;
                for (self.kernel_block[0..self.kernel_len]) |*v| v.* *= scale;
                @memset(self.kernel_block[self.kernel_len .. block_len * 2], 0);
                ffto.forward(self.kernel_block, work);
            }
        }
    };
}

pub fn BlockConvolver(comptime T: type) type {
    return struct {
        const Self = @This();
        const RFFT = fft.RealFFT(T);
        
        filter: *FIRFilter(T),
        fftin: *RFFT,
        ffto2: ?*RFFT,
        fftout: *RFFT,
        
        up_factor: usize,
        down_factor: usize,
        block_len2: usize,
        out_offset: usize,
        prev_input_len: usize,
        input_len: usize,
        latency_frac: T,
        latency: usize,
        
        up_shift: isize,
        down_shift: isize,
        input_delay: usize,
        
        prev_input: []T,
        cur_input: []T,
        cur_output: []T,
        
        in_data_left: usize,
        latency_left: usize,
        up_skip: usize,
        down_skip: usize,
        down_skip_init: usize,
        
        work_blocks: []T,
        /// FFT scratch (2 * block_len2), per convolver so the cached FFTs stay read-only.
        fft_work: []align(64) T,
        do_consume_latency: bool,
        cache: *FIRFilterCache,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator, filter: *FIRFilter(T), up_factor: usize, down_factor: usize, prev_latency: T, do_consume_latency: bool, cache: *FIRFilterCache) !*Self {
            std.debug.assert(up_factor > 0);
            std.debug.assert(down_factor > 0);
            std.debug.assert(prev_latency >= 0.0);

            // Owns the caller's reference to filter from here on, also on failure.
            errdefer cache.unref(filter);
            const self = try allocator.create(Self);
            errdefer allocator.destroy(self);
            
            self.* = .{
                .filter = filter,
                .up_factor = up_factor,
                .down_factor = down_factor,
                .do_consume_latency = do_consume_latency,
                .block_len2 = @as(usize, 2) << @as(std.math.Log2Int(usize), @intCast(filter.block_len_bits)),
                .work_blocks = &[_]T{},
                .fft_work = &.{},
                .prev_input = &[_]T{},
                .cur_input = &[_]T{},
                .cur_output = &[_]T{},
                .latency_frac = 0,
                .latency = 0,
                .out_offset = 0,
                .prev_input_len = 0,
                .input_len = 0,
                .up_shift = 0,
                .down_shift = 0,
                .input_delay = 0,
                .in_data_left = 0,
                .latency_left = 0,
                .up_skip = 0,
                .down_skip = 0,
                .down_skip_init = 0,
                .cache = cache,
                .allocator = allocator,
                .fftin = undefined,
                .ffto2 = null,
                .fftout = undefined,
            };
            
            var fftin_bits: usize = 0;
            self.up_shift = @as(isize, @intCast(base.getBitOccupancy(up_factor))) - 1;
            
            if ((@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(self.up_shift))) == up_factor) {
                fftin_bits = filter.block_len_bits + 1 - @as(usize, @intCast(self.up_shift));
                self.prev_input_len = (filter.kernel_len - 1 + up_factor - 1) / up_factor;
                self.input_len = self.block_len2 - self.prev_input_len * up_factor;
            } else {
                self.up_shift = -1;
                fftin_bits = filter.block_len_bits + 1;
                self.prev_input_len = filter.kernel_len - 1;
                self.input_len = self.block_len2 - self.prev_input_len;
            }

            self.out_offset = if (filter.is_zero_phase) filter.latency else 0;
            self.latency_frac = filter.latency_frac + prev_latency * @as(T, @floatFromInt(up_factor));
            self.latency = @as(usize, @intFromFloat(self.latency_frac));
            
            // Intentionally cast as isize for negative checks later
            const in_latency = @as(isize, @intCast(self.latency)) + @as(isize, @intCast(filter.latency)) - @as(isize, @intCast(self.out_offset));
            self.latency_frac -= @as(T, @floatFromInt(self.latency));
            self.latency_frac /= @as(T, @floatFromInt(down_factor));
            
            self.latency += self.input_len + filter.latency;

            var fftout_bits: usize = 0;
            self.input_delay = 0;
            self.down_skip_init = 0;
            self.down_shift = @as(isize, @intCast(base.getBitOccupancy(down_factor))) - 1;

            if ((@as(usize, 1) << @as(std.math.Log2Int(usize), @intCast(self.down_shift))) == down_factor) {
                fftout_bits = filter.block_len_bits + 1 - @as(usize, @intCast(self.down_shift));
                if (down_factor > 1) {
                    if (self.up_shift > 0) {
                        std.debug.assert(self.up_shift == 0);
                    } else {
                        const ilc = self.input_len & (down_factor - 1);
                        self.prev_input_len += ilc;
                        self.input_len -= ilc;
                        self.latency -= ilc;

                        const lc = in_latency & @as(isize, @intCast(down_factor - 1));
                        if (lc > 0) {
                            self.input_delay = down_factor - @as(usize, @intCast(lc));
                        }
                        if (!self.do_consume_latency) {
                            self.latency /= down_factor;
                        }
                    }
                }
            } else {
                fftout_bits = filter.block_len_bits + 1;
                self.down_shift = -1;
                
                if (!self.do_consume_latency and down_factor > 1) {
                    self.down_skip_init = self.latency % down_factor;
                    self.latency /= down_factor;
                }
            }

            self.fftin = try cache.getFFT(fftin_bits);

            if (fftout_bits == fftin_bits) {
                self.fftout = self.fftin;
            } else {
                const fftout = try cache.getFFT(fftout_bits);
                self.ffto2 = fftout;
                self.fftout = fftout;
            }

            // [cur_input | cur_output | fft_work | prev_input]: every part but the last stays 64-aligned.
            self.work_blocks = try allocator.alignedAlloc(T, .@"64", self.block_len2 * 4 + self.prev_input_len);
            @memset(self.work_blocks, 0);

            self.cur_input = self.work_blocks[0 .. self.block_len2];
            self.cur_output = self.work_blocks[self.block_len2 .. self.block_len2 * 2];
            self.fft_work = @alignCast(self.work_blocks[self.block_len2 * 2 .. self.block_len2 * 4]);
            self.prev_input = self.work_blocks[self.block_len2 * 4 ..];
            
            self.clear();
            return self;
        }

        pub fn deinit(self: *Self) void {
            // work_blocks is allocated 64-aligned (alignedAlloc above); the
            // field type is []T, so free it back through a matching aligned
            // slice or the allocator sees a 64-alloc / 8-free mismatch.
            if (self.work_blocks.len > 0) {
                const wb: []align(64) T = @alignCast(self.work_blocks);
                self.allocator.free(wb);
            }
            // fftin/ffto2 owned by FIRFilterCache.fft_cache — do not deinit here
            self.cache.unref(self.filter);
            self.allocator.destroy(self);
        }

        pub fn getInLenBeforeOutPos(self: *Self, req_out_pos: usize) usize {
            const upf = self.up_factor;
            const df = self.down_factor;
            return @intFromFloat(
                (@as(f64, @floatFromInt(self.latency)) + @as(f64, @floatFromInt(req_out_pos)) * @as(f64, @floatFromInt(df))) /
                    @as(f64, @floatFromInt(upf)) +
                    @as(f64, @floatCast(self.latency_frac)) * @as(f64, @floatFromInt(df)) / @as(f64, @floatFromInt(upf)),
            );
        }

        pub fn getLatency(self: *Self) usize {
            return if (self.do_consume_latency) 0 else self.latency;
        }

        pub fn getLatencyFrac(self: *Self) f64 {
            return @as(f64, @floatCast(self.latency_frac));
        }

        pub fn getMaxOutLen(self: *Self, max_in_len: usize) usize {
            return (max_in_len * self.up_factor + self.down_factor - 1) / self.down_factor + 1;
        }

        pub fn getName(ptr: *anyopaque) []const u8 {
            _ = ptr;
            return "BlockConvolver";
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

        fn processGeneric(ptr: *anyopaque, ip: [*]const f64, l: usize, op0: *[*]f64) usize {
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
            @memset(self.prev_input, 0);
            
            if (self.do_consume_latency) {
                self.latency_left = self.latency;
            } else {
                self.latency_left = 0;
                
                if (self.down_shift > 0) {
                    const shift = @as(std.math.Log2Int(usize), @intCast(self.down_shift));
                    @memset(self.cur_output[0 .. self.block_len2 >> shift], 0);
                } else {
                    @memset(self.cur_output[self.block_len2 - self.out_offset .. self.block_len2], 0);
                    @memset(self.cur_output[0 .. self.input_len - self.out_offset], 0);
                }
            }
            
            @memset(self.cur_input[0..self.input_delay], 0);
            
            self.in_data_left = self.input_len - self.input_delay;
            self.up_skip = 0;
            self.down_skip = self.down_skip_init;
        }

        fn copyUpsample(self: *Self, ip: *[]const T, op: []T, l0: usize) void {
            var input = ip.*;
            var b = @min(self.up_skip, l0);
            var remaining_l0 = l0;
            var op_idx: usize = 0;
            var ip_idx: usize = 0;

            if (b != 0) {
                self.up_skip -= b;
                remaining_l0 -= b;
                op[op_idx] = 0.0;
                op_idx += 1;
                while (b > 1) {
                    b -= 1;
                    op[op_idx] = 0.0;
                    op_idx += 1;
                }
            }

            const upf = self.up_factor;
            var l = remaining_l0 / upf;
            var lz = remaining_l0 - l * upf;

            if (upf == 3) {
                while (l != 0) {
                    op[op_idx + 0] = input[ip_idx];
                    op[op_idx + 1] = 0.0;
                    op[op_idx + 2] = 0.0;
                    ip_idx += 1;
                    op_idx += upf;
                    l -= 1;
                }
            } else if (upf == 5) {
                while (l != 0) {
                    op[op_idx + 0] = input[ip_idx];
                    op[op_idx + 1] = 0.0;
                    op[op_idx + 2] = 0.0;
                    op[op_idx + 3] = 0.0;
                    op[op_idx + 4] = 0.0;
                    ip_idx += 1;
                    op_idx += upf;
                    l -= 1;
                }
            } else {
                while (l != 0) {
                    op[op_idx] = input[ip_idx];
                    ip_idx += 1;
                    @memset(op[op_idx + 1 .. op_idx + upf], 0);
                    op_idx += upf;
                    l -= 1;
                }
            }

            if (lz != 0) {
                op[op_idx] = input[ip_idx];
                ip_idx += 1;
                op_idx += 1;
                self.up_skip = upf - lz;
                while (lz > 1) {
                    lz -= 1;
                    op[op_idx] = 0.0;
                    op_idx += 1;
                }
            }
            ip.* = input[ip_idx..];
        }

        fn copyToOutput(self: *Self, offs0: isize, output_ptr: *[]T, b0: usize, l0_ptr: *usize) void {
            var offs = offs0;
            var b = b0;

            if (offs < 0) {
                if (offs + @as(isize, @intCast(b)) <= 0) {
                    offs += @as(isize, @intCast(self.block_len2));
                } else {
                    self.copyToOutput(offs + @as(isize, @intCast(self.block_len2)), output_ptr, @as(usize, @intCast(-offs)), l0_ptr);
                    b = @as(usize, @intCast(@as(isize, @intCast(b)) + offs));
                    offs = 0;
                }
            }

            if (self.latency_left != 0) {
                if (self.latency_left >= b) {
                    self.latency_left -= b;
                    return;
                }
                offs += @as(isize, @intCast(self.latency_left));
                b -= self.latency_left;
                self.latency_left = 0;
            }

            const df = self.down_factor;
            var op0 = output_ptr.*;

            if (self.down_shift > 0) {
                var skip = @as(usize, @intCast(offs)) & (df - 1);
                if (skip > 0) {
                    skip = df - skip;
                    b -= skip;
                    offs += @as(isize, @intCast(skip));
                }
                if (b > 0) {
                    b = (b + df - 1) >> @as(std.math.Log2Int(usize), @intCast(self.down_shift));
                    const offs_u = @as(usize, @intCast(offs));
                    @memcpy(op0[0..b], self.cur_output[offs_u >> @as(std.math.Log2Int(usize), @intCast(self.down_shift)) .. (offs_u >> @as(std.math.Log2Int(usize), @intCast(self.down_shift))) + b]);
                    output_ptr.* = op0[b..];
                    l0_ptr.* += b;
                }
            } else {
                if (df > 1) {
                    const offs_u = @as(usize, @intCast(offs));
                    var ip_idx = offs_u + self.down_skip;
                    var l = (b + df - 1 - self.down_skip) / df;
                    self.down_skip = l * df + self.down_skip - b;
                    
                    l0_ptr.* += l;
                    var op_idx: usize = 0;
                    while (l > 0) {
                        op0[op_idx] = self.cur_output[ip_idx];
                        ip_idx += df;
                        op_idx += 1;
                        l -= 1;
                    }
                    output_ptr.* = op0[op_idx..];
                } else {
                    const offs_u = @as(usize, @intCast(offs));
                    @memcpy(op0[0..b], self.cur_output[offs_u .. offs_u + b]);
                    output_ptr.* = op0[b..];
                    l0_ptr.* += b;
                }
            }
        }

        fn mirrorInputSpectrum(self: *Self, p: []T) void {
            const bl1 = self.block_len2 >> @as(std.math.Log2Int(usize), @intCast(self.up_shift));
            const bl2 = bl1 + bl1;
            
            var i: usize = bl1 + 2;
            while (i < bl2) : (i += 2) {
                p[i] = p[bl2 - i];
                p[i + 1] = -p[bl2 - i + 1];
            }
            p[bl1] = p[1];
            p[bl1 + 1] = 0.0;
            p[1] = p[0];
            
            i = 1;
            while (i < @as(usize, @intCast(self.up_shift))) : (i += 1) {
                const z = bl1 << @as(std.math.Log2Int(usize), @intCast(i));
                @memcpy(p[z .. z + z], p[0..z]);
                p[z + 1] = 0.0;
            }
        }

        pub fn process(self: *Self, ip_ptr: []const T, op_full: []T) usize {
            var l = ip_ptr.len * self.up_factor;
            var ip = ip_ptr;
            var op = op_full;
            var l0: usize = 0;

            while (l > 0) {
                const offs = self.input_len - self.in_data_left;
                if (l < self.in_data_left) {
                    self.in_data_left -= l;
                    if (self.up_shift >= 0) {
                        const shift = @as(std.math.Log2Int(usize), @intCast(self.up_shift));
                        @memcpy(self.cur_input[offs >> shift .. (offs >> shift) + (l >> shift)], ip[0 .. l >> shift]);
                    } else {
                        self.copyUpsample(&ip, self.cur_input[offs..], l);
                    }
                    self.copyToOutput(@as(isize, @intCast(offs)) - @as(isize, @intCast(self.out_offset)), &op, l, &l0);
                    break;
                }
                
                const b = self.in_data_left;
                l -= b;
                self.in_data_left = self.input_len;
                var ilu: usize = undefined;
                
                if (self.up_shift >= 0) {
                    const shift = @as(std.math.Log2Int(usize), @intCast(self.up_shift));
                    const bu = b >> shift;
                    @memcpy(self.cur_input[offs >> shift .. (offs >> shift) + bu], ip[0..bu]);
                    ip = ip[bu..];
                    ilu = self.input_len >> shift;
                } else {
                    self.copyUpsample(&ip, self.cur_input[offs..], b);
                    ilu = self.input_len;
                }
                
                @memcpy(self.cur_input[ilu .. ilu + self.prev_input_len], self.prev_input[0..self.prev_input_len]);
                @memcpy(self.prev_input[0..self.prev_input_len], self.cur_input[ilu - self.prev_input_len .. (ilu - self.prev_input_len) + self.prev_input_len]);
                
                self.fftin.forward(self.cur_input, self.fft_work);
                if (self.up_shift > 0) {
                    self.mirrorInputSpectrum(self.cur_input);
                }
                
                if (self.filter.is_zero_phase) {
                    self.fftout.multiplyBlocksZP(self.filter.kernel_block, self.cur_input);
                } else {
                    self.fftout.multiplyBlocks(self.filter.kernel_block, self.cur_input);
                }
                
                if (self.down_shift > 0) {
                    const shift = @as(std.math.Log2Int(usize), @intCast(self.down_shift));
                    const z = self.block_len2 >> shift;
                    const kb = self.filter.kernel_block;
                    self.cur_input[1] = kb[z] * self.cur_input[z] - kb[z + 1] * self.cur_input[z + 1];
                }
                
                self.fftout.inverse(self.cur_input, self.fft_work);
                self.copyToOutput(@as(isize, @intCast(offs)) - @as(isize, @intCast(self.out_offset)), &op, b, &l0);
                
                const tmp = self.cur_input;
                self.cur_input = self.cur_output;
                self.cur_output = tmp;
            }
            return l0;
        }
    };
}

/// Filters, filter banks and FFT setups shared by every Resampler built on this cache.
/// With `initThreadSafe`, Resamplers on different threads may share it: init/deinit of
/// convolvers and interpolators take `mutex`; everything handed out is read-only after
/// publish, so process() never touches the cache. deinit only after all users are gone.
pub const FIRFilterCache = struct {
    const Self = @This();
    const FracDelayFilterBank = @import("interpolator.zig").FracDelayFilterBank;
    const max_cached: usize = 96;

    const Entry = struct {
        filter: *FIRFilter(f64),
        ref_count: usize,
        norm_freq: f64,
        trans_band: f64,
        atten: f64,
        phase: FilterPhaseResponse,
        gain: f64,
        next: ?*Entry,
    };

    allocator: std.mem.Allocator,
    head: ?*Entry,
    count: usize,
    fft_cache: fft.RealFFTCache,
    bank_cache: @import("interpolator.zig").FracDelayFilterBankCache,
    /// null: single-threaded, no locking.
    io: ?std.Io = null,
    mutex: std.Io.Mutex = .init,

    pub fn init(allocator: std.mem.Allocator) Self {
        return .{
            .allocator = allocator,
            .head = null,
            .count = 0,
            .fft_cache = fft.RealFFTCache.init(allocator),
            .bank_cache = @import("interpolator.zig").FracDelayFilterBankCache.init(allocator),
        };
    }

    /// Same as init, but every cache access locks one mutex through io.
    pub fn initThreadSafe(allocator: std.mem.Allocator, io: std.Io) Self {
        var self = init(allocator);
        self.io = io;
        return self;
    }

    fn lock(self: *Self) void {
        if (self.io) |io| self.mutex.lockUncancelable(io);
    }

    fn unlock(self: *Self) void {
        if (self.io) |io| self.mutex.unlock(io);
    }

    pub fn getFFT(self: *Self, len_bits: usize) !*fft.RealFFT(f64) {
        self.lock();
        defer self.unlock();
        return self.fft_cache.get(len_bits);
    }

    pub fn getBank(self: *Self, filter_fracs: i32, element_size: usize, interp_points: usize, req_atten: f64, is_third: bool, correction_gain: f64) !*FracDelayFilterBank {
        self.lock();
        defer self.unlock();
        return self.bank_cache.get(filter_fracs, element_size, interp_points, req_atten, is_third, correction_gain);
    }

    pub fn unrefBank(self: *Self, bank: *FracDelayFilterBank) void {
        self.lock();
        defer self.unlock();
        self.bank_cache.unref(bank);
    }

    pub fn deinit(self: *Self) void {
        var cur = self.head;
        while (cur) |entry| {
            const next = entry.next;
            entry.filter.deinit();
            self.allocator.destroy(entry);
            cur = next;
        }
        self.head = null;
        self.count = 0;
        self.fft_cache.deinit();
        self.bank_cache.deinit();
    }

    pub fn getLPFilter(self: *Self, norm_freq: f64, trans_band: f64, atten: f64, phase: FilterPhaseResponse, gain: f64) !*FIRFilter(f64) {
        self.lock();
        defer self.unlock();
        var prev: ?*Entry = null;
        var cur = self.head;

        while (cur) |entry| {
            if (entry.norm_freq == norm_freq and
                entry.trans_band == trans_band and
                entry.atten == atten and
                entry.phase == phase and
                entry.gain == gain)
            {
                entry.ref_count += 1;
                if (prev) |p| {
                    p.next = entry.next;
                    entry.next = self.head;
                    self.head = entry;
                }
                return entry.filter;
            }

            if (entry.next == null and self.count >= max_cached) {
                if (entry.ref_count == 0) {
                    if (prev) |p| { p.next = null; } else { self.head = null; }
                    entry.filter.deinit();
                    self.allocator.destroy(entry);
                    self.count -= 1;
                } else {
                    if (prev) |p| {
                        p.next = null;
                        entry.next = self.head;
                        self.head = entry;
                    }
                }
                break;
            }

            prev = entry;
            cur = entry.next;
        }

        const filter = try FIRFilter(f64).init(self.allocator, norm_freq, trans_band, atten, phase, gain);
        errdefer filter.deinit();
        const entry = try self.allocator.create(Entry);
        entry.* = .{
            .filter = filter,
            .ref_count = 1,
            .norm_freq = norm_freq,
            .trans_band = trans_band,
            .atten = atten,
            .phase = phase,
            .gain = gain,
            .next = self.head,
        };
        self.head = entry;
        self.count += 1;
        return filter;
    }

    pub fn unref(self: *Self, filter: *FIRFilter(f64)) void {
        self.lock();
        defer self.unlock();
        var cur = self.head;
        while (cur) |entry| {
            if (entry.filter == filter) {
                if (entry.ref_count > 0) entry.ref_count -= 1;
                return;
            }
            cur = entry.next;
        }
    }
};
