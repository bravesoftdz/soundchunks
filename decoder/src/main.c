#include <genesis.h>
#include <kdebug.h>

#include "../res/rsc.h"
#include "rsc.h"

#define RSC_SAMPLES_PER_CHUNK_SHIFT 2
#define RSC_SAMPLES_PER_CHUNK (1 << RSC_SAMPLES_PER_CHUNK_SHIFT)
#define RSC_CHUNKS_PER_FRAME 256

#define ROL8(x, n) ((x>>(8-n)) | (x<<n))

//#define USE_RSC_REF_DECODER
#define USE_YM_DEBUG_BIT

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
		*ymD0 = s;
	}

	inline auto void ymWT(void)
	{
#ifdef USE_YM_DEBUG_BIT
		while(!(*ymA0 & 0x40));
#else		
		while(!(*ymA0 & 0x01));
#endif		
	}
	
	inline auto void ymRT(void)
	{
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
	
	u8 track = 2;
	
	for (;;)
	{
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

		do
		{
			JOY_update();
		} while(JOY_readJoypad(0) & BUTTON_A);	

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

		// DAC init sample
		ymW0(0x2a, 0x80);
		ymW0(0x2c, 0x08);

#ifdef USE_YM_DEBUG_BIT
		// test
		ymW0(0x21, 0x41);
#else
		// Timer A enabled @ 26390Hz
		ymW0(0x24, 0xff);
		ymW0(0x25, 0x02);
		ymW0(0x27, 0x15);
#endif		

		u8 *chunks = &rsc[2];
		u16 blkCnt = 0, curChunkOff = 0, phase = 0;
		u8 sample, sampleLo, bitShift = 0, bsCtr = 0;
		vu8 vLo = 0xff;
		for (;;)
		{
			if (phase >= RSC_SAMPLES_PER_CHUNK * 0x100)
			{
				bitShift <<= 1;

				if (!bsCtr)
				{
					if (!blkCnt)
					{
						blkCnt = *rsc++;
						blkCnt |= (*rsc++) << 8;

						JOY_update();
						if (!blkCnt || (JOY_readJoypad(0) & BUTTON_A))
						{
							track++;
							break;
						}

						chunks = rsc;
						rsc += RSC_SAMPLES_PER_CHUNK * RSC_CHUNKS_PER_FRAME;
					}

					bitShift = *rsc++;
					bsCtr = 8;
					blkCnt--;
				}

				curChunkOff = *rsc++;
				phase = 0;
				bsCtr--;
			}

			sample = chunks[curChunkOff + phase];
			if (bitShift & 0x80)
			{
				sampleLo = (vLo << 3) & 8;
			}
			else
			{
				sampleLo = (sample << 3) & 8;
				sample = (sample >> 1) + 64;
			}
			
			ymWT();
			ymWSH(sample);
			ymWSL(sampleLo);
#ifndef USE_YM_DEBUG_BIT
			ymRT();
#endif			
			
			phase += 0x100;
		}
#else
		RSC_Init(rsc);

		for(;;)
		{
			if (JOY_readJoypad(0) & BUTTON_A)
			{
				VDP_drawText("Skipping...", 2, 14);
				RSC_StopTrack(FALSE);
			}

			if (RSC_IsTrackFinished())
			{
				RSC_Close();
				track++;
				break;
			}
			
			VDP_waitVSync();

			RSC_Set68kBusLockedFlag(TRUE);
			DMA_doDma(DMA_VRAM, 0, TILE_USER, 8 * 1024, 2);
			RSC_Set68kBusLockedFlag(FALSE);
		};
#endif
	}

	return 0;
}
