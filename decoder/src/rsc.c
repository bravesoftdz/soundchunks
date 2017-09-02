#include <genesis.h>

#include "rsc.h"
#include "z80_rsc.h"

#define Z80_RSCFLAGS 0x110
#define Z80_RSCBANK 0x111
#define Z80_68KLOCKED 0x112

static _voidCallback * vIntSysHandler;

static void RSC_vInt(void)
{
	RSC_Set68kBusLockedFlag(TRUE);
	vIntSysHandler();
	RSC_Set68kBusLockedFlag(FALSE);
}

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
	
	// install int handler to Set68kBusLockedFlag while 68k vInt is executing
	vIntSysHandler = internalVIntCB;
	internalVIntCB = RSC_vInt;
}

void RSC_Close(void)
{
	internalVIntCB = vIntSysHandler;	
	Z80_loadDriver(Z80_DRIVER_NULL, TRUE);
}

s8 RSC_IsTrackFinished(void)
{
	SYS_disableInts();
	Z80_requestBus(TRUE);
	u8 sts = Z80_read(Z80_RSCFLAGS);
	Z80_releaseBus();
	SYS_enableInts();
	
	return !!(sts & 4);
}

void RSC_Set68kBusLockedFlag(s8 flag)
{
	SYS_disableInts();
	flag = flag ? 0x80 : 0x00;	
	Z80_requestBus(TRUE);
	Z80_write(Z80_68KLOCKED, flag);
	Z80_releaseBus();
	SYS_enableInts();
}

void RSC_StopTrack(s8 waitFinished)
{
	// send stop request to Z80
	SYS_disableInts();
	Z80_requestBus(TRUE);
	Z80_write(Z80_RSCFLAGS, Z80_read(Z80_RSCFLAGS) | 8);
	Z80_releaseBus();
	SYS_enableInts();

	if (waitFinished)
	{
		do
		{
			waitSubTick(SUBTICKPERSECOND / 100);
		}
		while(!RSC_IsTrackFinished());
	}
}