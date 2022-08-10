# GliGli's SoundChunks audio codec

The goal of this project is to provide a novel lossy compressor of audio streams so that the decoder is so trivial to compute that even eg. MC68000 platforms can decode it with low CPU overhead.

So far, the encoder works reliably and a reference decoder is provided for further implementation help.

At full bitrate, quality is slightly better than ADPCM, and bitrate is lower.

Project page: https://github.com/gligli/soundchunks (source code, releases)

_Author: GliGli_ / _License: GNU GPL v3_

## File Format

The .gsc **stream** is one or more **frames** with a header each. Inside are **chunks** (which are short samples of the same size and have an integer attenuation), and **chunk indexes** (with attributes like 'negative' and 'reversed').
The decoder outputs sound by parsing the **chunk indexes** list of a **frame**, and outputing corresponding attenuated **chunks** corrected with the attributes.

### Frames
All multi-byte elements are **Little-Endian**
- StreamVersion (1 byte) : **Stream** binary format version.
- ChannelCount (1 byte) : Audio channel count (eg: 2 for stereo).
- ChunkCount (2 bytes) : Number of **chunks** of ChunkSize sole samples in the **frame** (eg: 4096 at full bitrate).
- ChunkBitDepth (1 byte) : **Chunks** format; 8: Unsigned 8 bits sole samples; 12: Unsigned 12 bits sole samples.
- ChunkSize (1 byte) : Number of sole samples in a **chunk**;
- SampleRate (3 bytes) : Audio sample rate (eg: 44100).
- ChunkBlend (1 byte) : Reserved (should be 0).
- AttenuationDivider (2 bytes) : 1 / AttenuationDivider is the increment for **chunks** attenuations (eg: chunk attenuation = 2 -> attenuation = 1 / AttenuationDivider + 2 / AttenuationDivider + 3 / AttenuationDivider).
- Attenuations (4 bits * ChunkCount; byte padded) : Attenuations, one per **chunk**. This is how many times 'Increment / AttenuationDivider' should be accumulated and multiplied to each sole sample to give an usable **chunk**.
- Chunk sole samples ((ChunkBitDepth bits per sole sample) * ChunkSize * ChunkCount elements; word padded for 12 bits): Actual *chunks* of ChunkSize sole samples.
- FrameLength (4 bytes) : Number of **chunk indexes** in the **frame**.
- ChunkIndexes ((3-4 bits + variable length coded) * FrameLength; word padded) : 
  * Negative (1 bit) : The **chunk** should be played with all sole samples in two's complement or with negative phase.
  * Reversed (1 bit; only if StreamVersion > 0) : The **chunk** should be played with all sole samples in backwards order.
  * VariableCodingHeader (2 bits) : 0b00: chunk index is 3 bits; 0b01: chunk index is 6 bits; 0b10: chunk index is 9 bits; 0b11: chunk index is 12 bits; 
  * Chunk index (variable length coded; cf. VariableCodingHeader): *Chunk index* from the frame.
