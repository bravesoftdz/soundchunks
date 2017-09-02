// MegaDrive 68k API for replaying RSE SoundChunks PCM streams with the Z80
// By GliGli

#pragma once

void RSC_Init(const u8 * rsc_track);
void RSC_Close(void);
void RSC_Set68kBusLockedFlag(s8 flag);
s8 RSC_IsTrackFinished(void);
void RSC_StopTrack(s8 waitFinished);
