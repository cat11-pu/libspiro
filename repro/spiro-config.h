/* Minimal config for standalone repro build (no autotools needed). */
#ifndef SPIRO_CONFIG_MIN_H
#define SPIRO_CONFIG_MIN_H
#define HAVE_ISFINITE 1
#define LS_VERSION_STR "repro-build"
/* Intentionally do NOT define CHECK_INPUT_FINITENESS: default build. */
/* Mirrors configure.ac AH_BOTTOM definition. */
#if HAVE_ISFINITE
#define IS_FINITE(x) isfinite(x)
#else
#if HAVE_FINITE
#define IS_FINITE(x) finite(x)
#endif
#endif
#endif
