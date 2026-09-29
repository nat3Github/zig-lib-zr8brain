const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zpffft = b.dependency("zpffft", .{ .target = target, .optimize = optimize });

    const mod = b.addModule("r8brain", .{
        .root_source_file = b.path("zig-src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    mod.addImport("zpffft", zpffft.module("pffft"));
}
