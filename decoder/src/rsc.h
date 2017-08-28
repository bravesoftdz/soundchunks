// MegaDrive 68k API for replaying RSE SoundChunks PCM streams with the Z80
// By GliGli

#pragma once

void RSC_Init(const u8 * rsc_track);
s8 RSC_IsTrackFinished(void);
void RSC_StopTrack(s8 waitFinished);
