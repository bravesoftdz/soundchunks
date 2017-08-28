#include <genesis.h>
#include <kdebug.h>

int main()
{
	//Z80_startReset();
	Z80_requestBus(TRUE);
	
	vu8 *ymA0 = (vu8 *) 0xa04000;
	vu8 *ymD0 = (vu8 *) 0xa04001;
	vu8 *ymA1 = (vu8 *) 0xa04002;
	vu8 *ymD1 = (vu8 *) 0xa04003;
	
	inline auto void ymW0(u8 a, u8 d)
	{
		if (a)
		{
			*ymA0 = a;
		}
		*ymD0 = d;
		while(*ymA0 & 0x80);
	}
	
	inline auto void ymWS(u8 s)
	{
		*ymA0 = 0x2a;
		*ymD0 = s;
	}

	inline auto void ymWT(void)
	{
		while(!(*ymA0 & 0x01));
		*ymA0 = 0x27;
		*ymD0 = 0x15;
	}
	
	inline auto void ymW1(u8 a, u8 d)
	{
		if (a)
		{
			*ymA1 = a;
		}
		*ymD1 = d;
		while(*ymA0 & 0x80);
	}
	
	// Enable DAC
	ymW0(0x2b, 0x80);
	
	// DAC Pan = LR
	ymW1(0xb6, 0xc0);

	// Timer A enabled @ 26390Hz
	ymW0(0x24, 0xff);
	ymW0(0x25, 0x02);
	ymW0(0x27, 0x15);

	// 9th DAC bit = 1
	ymW0(0x2c, 0x08);

	// DAC middle sample
	ymW0(0x2a, 0x80);
	
	u16 sine[256];
	for(s16 i = 0; i < 256 ; ++i)
	{
		fix32 f =min((sinFix32(i << 2) >> 3) + 128, 255);
		
//		if (f >= 128) f = max(128, f - 3);
		sine[i] = f;
	}
	
	VDP_drawText("init", 0, 0);
	
	SYS_disableInts();
	
	u8 ph = 0;
	for(;;)
	{
		ymWS(sine[ph++]);
		ymWT();
	}
	return 0;
}
