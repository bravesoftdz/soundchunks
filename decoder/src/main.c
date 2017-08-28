#include <genesis.h>
#include <kdebug.h>

#include "../res/rsc.h"

#define RSC_SAMPLES_PER_CHUNK_SHIFT 2
#define RSC_SAMPLES_PER_CHUNK (1 << RSC_SAMPLES_PER_CHUNK_SHIFT)
#define RSC_CHUNKS_PER_FRAME 256

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
	
	VDP_drawText("RSE SoundChunks decoder", 4, 4);
	
start:	
	SYS_disableInts();

	u8 * rsc = (u8 *) rsc_data;	
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
				if (!idxCnt) goto start;

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

#if 0
		KLog_U4("bsCtr ", bsCtr, " phase ", phase, " idxCnt ", idxCnt, " bitShift ", bitShift);
		VDP_waitVSync();
#endif		
	}
	return 0;
}
