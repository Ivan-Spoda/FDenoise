# FDenoise:

This is a batch audio noise removal tool, can reduce noise/artifacts by volume treshold.
Supported I/O formats: `wav`, `flac`, `aiff`.

**Examples**: [Original](https://github.com/Ivan-Spoda/FDenoise/tree/master/examples/Original) | [Processed](https://github.com/Ivan-Spoda/FDenoise/tree/master/examples/Processed). For velocity layers used `X_021_(velocity_layer)`.

### Settings:
| **Argument**                 | **Value**  |
|:-----------------------------|:-----------|
| **FFT Size**                 | `3`        |
| **Min Amplitude**            | `-90`      |
| **Max Amplitude**            | `-60`      |
| **Treshold**                 | `-90`      |
| **Format**                   | `flac`     |
| **Bit Depth**                | `24`       |
| **Resolution**               | `8`        |

### Usage:

Open CMD, go to FDenoise.exe path, run: `FDenoise.exe --help` for get addition instructions. \
**Example**: 
```bash
FDenoise.exe -i "inputPath" -o "outputPath" --fs 3 --ampMin -90 --ampMax -60 -t -90 --format flac --bit 24 --res 8
```

### Args table

| **Argument**   | **Description**                                                                        |
|:---------------|:---------------------------------------------------------------------------------------|
| `-i`           | Input path.                                                                            |
| `-o`           | Out path.                                                                              |
| `-t`           | Filtering bins volume treshold (dB).                                                   |
| `--fs`         | FFT Size (8192, 16384, 32768, 65536, 131072).                                          |
| `--res`        | FFT Resolution multiplier/Padding (1, 2, 4, 8).                                        |
| `--to`         | Time Overlap Divisor (2, 4, 8, 16, 32, 64, 128).                                       |
| `--ampMin`     | Make signal type contrast (dB).                                                        |
| `--ampMax`     | Make signal type contrast (dB).                                                        |
| `--format`     | Change output format (flac, wav, aiff).                                                |
| `--bit`        | Change output bit-depth (16, 24, 32).                                                  |

### Thirdparty libs

- [fftw3](https://fftw.org/)
- [libsndfile](https://libsndfile.github.io/libsndfile/)
