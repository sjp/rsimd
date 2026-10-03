/* Functions every tier translation unit defines, independent of kernels. */

/* The tier's id, as in rsimd_tier_names; confirms the object was linked. */
const char *RSIMD_KERNEL(tier_name)(void);
const char *RSIMD_KERNEL(tier_name)(void) { return RSIMD_TIER_STRING; }

#if RSIMD_TIER_IS(sve)
/* SVE vector length in bits. Only called when the CPU supports SVE. */
int rsimd_sve_vl_bits_sve(void);
int rsimd_sve_vl_bits_sve(void) { return (int) svcntb() * 8; }
#endif
