/* Functions every tier translation unit defines, independent of kernels. */

/* Kernel for the tier_name slot: the tier's id, as in rsimd_tier_names. */
const char *RSIMD_KERNEL(tier_name)(void);
const char *RSIMD_KERNEL(tier_name)(void) { return RSIMD_TIER_STRING; }

/* Kernel for the fill_probe slot, which only the none tier implements so
   that the dispatcher's fill-down is always in use. */
#if RSIMD_TIER_IS(none)
const char *RSIMD_KERNEL(fill_probe)(void);
const char *RSIMD_KERNEL(fill_probe)(void) { return RSIMD_TIER_STRING; }
#else
#define RSIMD_SKIP_fill_probe 1
#endif

#if RSIMD_TIER_IS(sve)
/* SVE vector length in bits. Only called when the CPU supports SVE. */
int rsimd_sve_vl_bits_sve(void);
int rsimd_sve_vl_bits_sve(void) { return (int) svcntb() * 8; }
#endif
