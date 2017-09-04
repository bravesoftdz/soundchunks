#include <genesis.h>
#include <kdebug.h>

#include "../res/rsc.h"
#include "rsc.h"

#define RSC_SAMPLES_PER_CHUNK_SHIFT 2
#define RSC_SAMPLES_PER_CHUNK (1 << RSC_SAMPLES_PER_CHUNK_SHIFT)
#define RSC_CHUNKS_PER_FRAME 256

#define ROL8(x, n) ((x>>(8-n)) | (x<<n))

//#define USE_RSC_REF_DECODER
//#define USE_YM_DEBUG_BIT

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
	
	u8 track = 0;
	
	for (;;)
	{
		VDP_clearPlan(PLAN_A, TRUE);
		VDP_drawText("RSE SoundChunks decoder", 8, 3);
#ifdef USE_RSC_REF_DECODER
		VDP_drawText("(68k ref mode)", 12, 4);
#else
		VDP_drawText("(Z80 mode)", 14, 4);
#endif	
		VDP_drawText("(Press A for next track)", 2, 22);

		u8 * rsc = NULL;

		track = track % 3;

		do
		{
			JOY_update();
		} while(JOY_readJoypad(0) & BUTTON_A);	

		switch(track)
		{
		case 0:
//			VDP_drawText("Kavinsky - Outrun Prelude", 2, 20);
			VDP_drawText("Queen - We Will Rock You", 2, 20);
			rsc = (u8 *) rsc_1;	
			break;
		case 1:
			VDP_drawText("David Whittaker - Shadow of the Beast", 2, 20);
			rsc = (u8 *) rsc_2;
			break;
		case 2:
//			VDP_drawText("Wintergatan - Live Intro", 2, 20);
			VDP_drawText("Robots with Rayguns - Body Motion", 2, 20);
			rsc = (u8 *) rsc_3;
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

		u8 *chunks = &rsc[1];
		u16 blkCnt = 0, curChunkOff = 0, phase = RSC_SAMPLES_PER_CHUNK * 0x100;
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
						blkCnt = (*rsc++) << 8;

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
		u16 tileScr = 0;
		s16 vscrolla[32];
		s16 vscrollb[32];
		
		auto void hInt(void)
		{
			static s16 avg = 0;
			avg = ((avg << 1) + RSC_GetLastSample()) >> 2;
			vscrollb[tileScr++] = avg;
		}
		
		RSC_Init(rsc);

		memsetU16((u16 *)vscrolla, 0, 32);
		memsetU16((u16 *)vscrollb, 0, 32);
		
		VDP_setScrollingMode(HSCROLL_PLANE, VSCROLL_2TILE);
		for (int y = 0; y < screenHeight / 8; ++y)
			for (int x = 0; x < screenWidth / 8; ++x)
				VDP_setTileMapXY(PLAN_B, TILE_ATTR_FULL(PAL3, FALSE, FALSE, FALSE, (y == 12) ? 15 : 5 - abs(12 - y)), x, y);
		
		internalHIntCB = hInt;
		VDP_setHIntCounter(10);
		VDP_setHInterrupt(TRUE);

		u8 pos = 0;
		for(;;)
		{
			if (JOY_readJoypad(0) & BUTTON_A)
			{
				VDP_drawText("Skipping...", 2, 12);
				RSC_StopTrack(FALSE);
			}

			if (RSC_IsTrackFinished())
			{
				RSC_Close();
				track++;
				break;
			}
			
			VDP_waitVSync();

			tileScr = 0;
			
			VDP_setEnable(FALSE);
			RSC_Set68kBusLockedFlag(TRUE);

			DMA_doDma(DMA_VRAM, 0, TILE_USER, (((pos++) & 0x3) + 1)  * 1024, 2);

			VDP_setVerticalScrollTile(PLAN_A, 0, vscrolla, 20, TRUE);
			VDP_setVerticalScrollTile(PLAN_B, 0, vscrollb, 20, TRUE);

			RSC_Set68kBusLockedFlag(FALSE);
			VDP_setEnable(TRUE);
			
		};
#endif
	}

	return 0;
}
