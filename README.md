# zr8brain

Pure-Zig port of [r8brain-free-src](https://github.com/avaneev/r8brain-free-src) — a high-quality professional audio sample rate converter by Aleksey Vaneev of Voxengo.

Converts audio between arbitrary sample rates (integer and fractional). Uses 2× oversampling + polynomial-interpolated sinc fractional delay filters. Power-of-2 and 3×power-of-2 ratios take a fast path automatically.

## Quick Start

```sh
zig fetch --save=r8brain git+https://github.com/nat3Github/zig-lib-zr8brain#<commit>
```

Import in `build.zig`:

```zig
const r8brain_dep = b.dependency("r8brain", .{ .optimize = optimize, .target = target });
module.addImport("r8brain", r8brain_dep.module("r8brain"));
```

Basic usage:

```zig
const r8b = @import("r8brain");

var cache = r8b.FIRFilterCache.init(allocator);
defer cache.deinit();

var resampler = try r8b.R8bResampler.init(
    allocator,
    44100.0,   // input rate
    48000.0,   // output rate
    1024,      // max input frames per call
    2.0,       // transition band (% of Nyquist)
    &cache,
);
defer resampler.deinit();

var out_buf: [2048]f32 = undefined;
const produced = resampler.process(input_f32, &out_buf);
```

## API

### `R8bResampler`

High-level f32 wrapper. Internally processes in f64.

| Function | Description |
|----------|-------------|
| `init(alloc, in_rate, out_rate, max_in_len, trans_band, cache)` | Create with `.quality24` (180 dB) |
| `initWithQuality(..., quality, cache)` | Create with explicit quality preset |
| `process(input: []const f32, output: []f32) usize` | Resample; returns samples written |
| `clear()` | Reset resampler state |
| `getLatency() i32` | Integer latency in output samples |
| `getLatencyFrac() f64` | Fractional latency |
| `getMaxOutLen(max_in_len: i32) i32` | Max output for given input length |
| `getInLenBeforeOutPos(req_out_pos: i32) i32` | Input needed before output position |
| `deinit()` | Free resources |

### `ResamplerQuality`

| Preset | Attenuation | Equivalent C++ class |
|--------|-------------|----------------------|
| `.quality16` | 136 dB | `CDSPResampler16` |
| `.quality16ir` | 110 dB | `CDSPResampler16IR` |
| `.quality24` | 180 dB | `CDSPResampler24` (default) |

### Shared Caches

Reuse these across multiple resamplers to avoid redundant filter computation:

| Type | Caches |
|------|--------|
| `FIRFilterCache` | FIR low-pass filters |
| `RealFFTCache` | FFT plan objects |
| `FracDelayFilterBankCache` | Fractional delay filter banks |

## Dependencies

[zig-lib-zpffft](https://github.com/nat3Github/zig-lib-zpffft) (pure-Zig PFFFT) — nothing else, no libc.
All memory comes from the allocator you pass in; no globals, no threads.

## Credits

Algorithm and original implementation by **Aleksey Vaneev of Voxengo**.
Per the original library's request: "Sample rate converter designed by Aleksey Vaneev of Voxengo".

## License

MIT (same as the original r8brain-free-src), see [LICENSE](LICENSE). The original C++ README is in [README_r8brain.md](README_r8brain.md).
