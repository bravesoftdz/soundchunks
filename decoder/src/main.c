#include <genesis.h>
#include <kdebug.h>

#include "../res/rsc.h"
#include "rsc.h"

#define RSC_SAMPLES_PER_CHUNK_SHIFT 2
#define RSC_SAMPLES_PER_CHUNK (1 << RSC_SAMPLES_PER_CHUNK_SHIFT)
#define RSC_CHUNKS_PER_FRAME 256

#define USE_RSC_REF_DECODER

int main()
{
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
	
	inline auto void ymWSH(u8 s)
	{
		*ymA0 = 0x2a;
		*ymD0 = s;
	}

	inline auto void ymWSL(u8 s)
	{
		*ymA0 = 0x2c;
		*ymD0 = s << 3;
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
	
	u8 track = 0;
	
start:	
	VDP_clearPlan(PLAN_A, TRUE);
	VDP_drawText("RSE SoundChunks decoder", 8, 4);
#ifdef USE_RSC_REF_DECODER
	VDP_drawText("(68k ref mode)", 12, 5);
#else
	VDP_drawText("(Z80 mode)", 14, 5);
#endif	
	VDP_drawText("(Press A for next track)", 2, 12);
	
	u8 * rsc = NULL;

	track = track % 3;
	
	switch(track)
	{
	case 0:
		VDP_drawText("Kavinsky - Outrun Prelude", 2, 10);
		rsc = (u8 *) rsc_kav;	
		break;
	case 1:
		VDP_drawText("David Whittaker - Shadow of the Beast", 2, 10);
		rsc = (u8 *) rsc_sob;
		break;
	case 2:
		VDP_drawText("Wintergatan - Live Intro", 2, 10);
		rsc = (u8 *) rsc_win;
		break;
	}
	
#ifdef USE_RSC_REF_DECODER
	SYS_disableInts();
	Z80_requestBus(TRUE);

	// Enable DAC
	ymW0(0x2b, 0x80);
	
	// DAC Pan = LR
	ymW1(0xb6, 0xc0);

	// Timer A enabled @ 26390Hz
	ymW0(0x24, 0xff);
	ymW0(0x25, 0x02);
	ymW0(0x27, 0x15);

	// DAC init sample
	ymW0(0x2c, 0x08);
	ymW0(0x2a, 0x80);
	
	do
	{
		JOY_update();
	} while(JOY_readJoypad(0) & BUTTON_A);	

	u8 *chunks = &rsc[2];
	u16 idxCnt = 0, bsCtr = 0, phase = 0, curChunkOff = 0;
	u8 bitShift, sample;
	for(;;)
	{
		if(!phase)
		{
			if (!idxCnt)
			{
				idxCnt = *rsc++;
				idxCnt |= (*rsc++) << 8;
				
				JOY_update();
				if (!idxCnt || (JOY_readJoypad(0) & BUTTON_A))
				{
					track++;
					goto start;
				}

				chunks = rsc;
				rsc += RSC_SAMPLES_PER_CHUNK * RSC_CHUNKS_PER_FRAME;
				bsCtr = 0;
				phase = 0;
			}
			
			bitShift <<= 1;

			if (!bsCtr)
			{
				bitShift = *rsc++;
				bsCtr = 8;
			}

			curChunkOff = (*rsc++) << RSC_SAMPLES_PER_CHUNK_SHIFT;
			phase = RSC_SAMPLES_PER_CHUNK;
			bsCtr--;
			idxCnt--;
		}
		
		phase--;
		
		sample = chunks[curChunkOff + phase];
		if (bitShift & 0x80)
		{
			ymWSH(sample);
			ymWSL(1);
		}
		else
		{
			ymWSH((sample >> 1) + 64);
			ymWSL(sample & 1);
		}

		ymWT();
	}
#else
	RSC_Init(rsc);
	
	for(;;)
	{
		if (JOY_readJoypad(0) & BUTTON_A)
		{
			RSC_StopTrack(FALSE);
		}
		
		if (RSC_IsTrackFinished())
		{
			track++;
			goto start;
		}
		
		VDP_waitVSync();
	}
#endif

	return 0;
}
