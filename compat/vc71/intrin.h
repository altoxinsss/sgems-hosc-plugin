#ifndef WITHSGEMS_VC71_INTRIN_COMPAT_H
#define WITHSGEMS_VC71_INTRIN_COMPAT_H

/*
 * The VC++ Toolkit 2003 compiler implements the x86 intrinsics below but
 * does not ship the later intrin.h declaration header.  The legacy SGeMS SDK
 * carries a modern Boost header which includes intrin.h unconditionally for
 * all Microsoft compilers.  Keep this declaration-only shim in the VC7.1
 * build path so the host-compatible compiler remains self-contained.
 */

#ifdef __cplusplus
extern "C" {
#endif

unsigned char __cdecl _BitScanForward(unsigned long* index, unsigned long mask);
unsigned char __cdecl _BitScanReverse(unsigned long* index, unsigned long mask);
unsigned long __cdecl _byteswap_ulong(unsigned long value);
unsigned __int64 __cdecl _byteswap_uint64(unsigned __int64 value);

#ifdef __cplusplus
}
#endif

#pragma intrinsic(_BitScanForward)
#pragma intrinsic(_BitScanReverse)
#pragma intrinsic(_byteswap_ulong)
#pragma intrinsic(_byteswap_uint64)

#endif
