#include <genesis.h>

#include "rsc.h"
#include "z80_rsc.h"

#define Z80_RSCFLAGS 0x110
#define Z80_RSCBANK 0x111
#define Z80_68KLOCKED 0x112

void RSC_Init(const u8 * rsc_track)
{
	u8 sts;
	
	Z80_loadCustomDriver(z80_rsc, sizeof(z80_rsc));

	Z80_requestBus(TRUE);
	
	// config
	Z80_write(Z80_RSCBANK, (u32)rsc_track >> 15);
	
	// 68k side init done
	Z80_write(Z80_RSCFLAGS, 0x01);
	Z80_releaseBus();

	// wait Z80 size init done
	do
	{
		waitSubTick(SUBTICKPERSECOND / 100);
		Z80_requestBus(TRUE);
		sts = Z80_read(Z80_RSCFLAGS);
		Z80_releaseBus();
	}
	while(!(sts & 2));
}

s8 RSC_IsTrackFinished(void)
{
	Z80_requestBus(TRUE);
	u8 sts = Z80_read(Z80_RSCFLAGS);
	Z80_releaseBus();
	
	return !!(sts & 4);
}

void RSC_Set68kBusLockedFlag(s8 flag)
{
	flag = flag ? 0x80 : 0x00;	
	Z80_requestBus(TRUE);
	Z80_write(Z80_68KLOCKED, flag);
	Z80_releaseBus();
}

void RSC_StopTrack(s8 waitFinished)
{
	// send stop request to Z80
	Z80_requestBus(TRUE);
	Z80_write(Z80_RSCFLAGS, Z80_read(Z80_RSCFLAGS) | 8);
	Z80_releaseBus();

	if (waitFinished)
	{
		do
		{
			waitSubTick(SUBTICKPERSECOND / 100);
		}
		while(!RSC_IsTrackFinished());
	}
}