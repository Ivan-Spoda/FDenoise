#include <QCommandLineParser>
#include <QCoreApplication>
#include <QDebug>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QMutex>
#include <QVector>
#include <QtConcurrent>
#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <sndfile.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static QMutex fftwPlanMutex;

struct ProcessSettings
{
    QString inputPath;
    QString outputPath;
    int windowSize;
    int resolutionFactor;
    int hopDivisor;
    double thresholdDB;
    int outputFormat;
    int outputBitDepth;
};

struct FileJob
{
    QString inputPath;
    ProcessSettings settings;
};

std::vector<double> makeHannWindow(int size)
{
    std::vector<double> window(size);
    for (int i = 0; i < size; ++i) {
        window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (size - 1)));
    }
    return window;
}

inline double linearToDb(double linear)
{
    if (linear < 1e-20)
        return -400.0;
    return 20.0 * std::log10(linear);
}

void processChannel(const std::vector<double> &rawInput,
                    std::vector<double> &output,
                    const ProcessSettings &settings)
{
    int winSize = settings.windowSize;
    int padding = winSize;

    std::vector<double> input(padding, 0.0);
    input.insert(input.end(), rawInput.begin(), rawInput.end());
    input.insert(input.end(), padding, 0.0);

    int fftSize = winSize * settings.resolutionFactor;
    int hopSize = winSize / settings.hopDivisor;

    if (hopSize <= 0)
        hopSize = winSize / 4;

    fftwPlanMutex.lock();
    double *inBuffer = (double *) fftw_malloc(sizeof(double) * fftSize);
    fftw_complex *outSpec = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (fftSize / 2 + 1));
    fftw_plan pFwd = fftw_plan_dft_r2c_1d(fftSize, inBuffer, outSpec, FFTW_ESTIMATE);
    fftw_plan pInv = fftw_plan_dft_c2r_1d(fftSize, outSpec, inBuffer, FFTW_ESTIMATE);
    fftwPlanMutex.unlock();

    std::vector<double> window = makeHannWindow(winSize);
    std::vector<double> tempOutput(input.size() + fftSize, 0.0);
    std::vector<double> weightBuffer(tempOutput.size(), 0.0);

    double fftScale = 1.0 / (double) fftSize;
    int numFrames = (input.size() - winSize) / hopSize;
    if (numFrames < 0)
        numFrames = 0;

    for (int i = 0; i <= numFrames; ++i) {
        int pos = i * hopSize;

        std::fill(inBuffer, inBuffer + fftSize, 0.0);

        for (int j = 0; j < winSize; ++j) {
            if (pos + j < input.size()) {
                inBuffer[j] = input[pos + j] * window[j];
            }
        }

        fftw_execute(pFwd);

        int bins = fftSize / 2 + 1;
        for (int k = 0; k < bins; ++k) {
            double re = outSpec[k][0];
            double im = outSpec[k][1];

            double mag = std::sqrt(re * re + im * im);
            double magNorm = mag * (2.0 / winSize);
            double db = linearToDb(magNorm);

            if (db < settings.thresholdDB) {
                outSpec[k][0] = 0.0;
                outSpec[k][1] = 0.0;
            }
        }

        fftw_execute(pInv);

        for (int j = 0; j < winSize; ++j) {
            if (pos + j < tempOutput.size()) {
                double sample = inBuffer[j] * fftScale;

                tempOutput[pos + j] += sample * window[j];
                weightBuffer[pos + j] += window[j] * window[j];
            }
        }
    }

    for (size_t i = 0; i < tempOutput.size(); ++i) {
        if (weightBuffer[i] > 1e-9) {
            tempOutput[i] /= weightBuffer[i];
        }
    }

    output.reserve(rawInput.size());
    for (size_t i = 0; i < rawInput.size(); ++i) {
        if (padding + i < tempOutput.size()) {
            output.push_back(tempOutput[padding + i]);
        } else {
            output.push_back(0.0);
        }
    }

    fftwPlanMutex.lock();
    fftw_destroy_plan(pFwd);
    fftw_destroy_plan(pInv);
    fftw_free(inBuffer);
    fftw_free(outSpec);
    fftwPlanMutex.unlock();
}

void processFileJob(const FileJob &job)
{
    QFileInfo fi(job.inputPath);
    const ProcessSettings &settings = job.settings;

    QString ext;
    int typeMask = settings.outputFormat & SF_FORMAT_TYPEMASK;
    if (typeMask == SF_FORMAT_FLAC)
        ext = ".flac";
    else if (typeMask == SF_FORMAT_AIFF)
        ext = ".aiff";
    else
        ext = ".wav";

    QString outName = settings.outputPath + QDir::separator() + fi.completeBaseName() + ext;

#ifdef _WIN32
    QByteArray inPathB = job.inputPath.toLocal8Bit();
    QByteArray outPathB = outName.toLocal8Bit();
#else
    QByteArray inPathB = job.inputPath.toUtf8();
    QByteArray outPathB = outName.toUtf8();
#endif

    SF_INFO sfInfoIn;
    std::memset(&sfInfoIn, 0, sizeof(sfInfoIn));
    SNDFILE *inFile = sf_open(inPathB.constData(), SFM_READ, &sfInfoIn);

    if (!inFile) {
        qCritical() << "Error opening:" << fi.fileName();
        return;
    }

    qInfo() << "Processing:" << fi.fileName() << "| Res: x" << settings.resolutionFactor
            << "| Win:" << settings.windowSize << "| Threshold:" << settings.thresholdDB << "dB";

    sf_count_t numFrames = sfInfoIn.frames;
    if (numFrames == 0) {
        sf_close(inFile);
        return;
    }

    std::vector<double> interleavedData(numFrames * sfInfoIn.channels);
    sf_count_t readFrames = sf_readf_double(inFile, interleavedData.data(), numFrames);
    sf_close(inFile);

    int channels = sfInfoIn.channels;
    std::vector<std::vector<double>> channelData(channels);
    for (int c = 0; c < channels; ++c) {
        channelData[c].reserve(readFrames);
    }

    for (size_t i = 0; i < (size_t) readFrames * channels; i += channels) {
        for (int c = 0; c < channels; ++c) {
            channelData[c].push_back(interleavedData[i + c]);
        }
    }

    std::vector<std::vector<double>> processedChannels(channels);

    for (int c = 0; c < channels; ++c) {
        processChannel(channelData[c], processedChannels[c], settings);
    }

    size_t maxLen = 0;
    for (int c = 0; c < channels; ++c) {
        if (processedChannels[c].size() > maxLen)
            maxLen = processedChannels[c].size();
    }

    std::vector<double> outInterleaved;
    outInterleaved.reserve(maxLen * channels);

    for (size_t i = 0; i < maxLen; ++i) {
        for (int c = 0; c < channels; ++c) {
            if (i < processedChannels[c].size())
                outInterleaved.push_back(processedChannels[c][i]);
            else
                outInterleaved.push_back(0.0);
        }
    }

    SF_INFO sfInfoOut;
    std::memset(&sfInfoOut, 0, sizeof(sfInfoOut));
    sfInfoOut.channels = channels;
    sfInfoOut.samplerate = sfInfoIn.samplerate;
    sfInfoOut.format = (settings.outputFormat & SF_FORMAT_TYPEMASK)
                       | (settings.outputBitDepth & SF_FORMAT_SUBMASK);

    if (!sf_format_check(&sfInfoOut)) {
        sfInfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;
    }

    SNDFILE *outFile = sf_open(outPathB.constData(), SFM_WRITE, &sfInfoOut);
    if (!outFile) {
        qCritical() << "Write error:" << outName;
        return;
    }

    sf_writef_double(outFile, outInterleaved.data(), maxLen);
    sf_close(outFile);

    qInfo() << "Saved:" << outName;
}

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName(TARGET_NAME);
    QCoreApplication::setApplicationVersion(TARGET_VER);

    QCommandLineParser parser;
    parser.addHelpOption();

    parser.addOption({{"i", "input"}, "Input path.", "path"});
    parser.addOption({{"o", "output"}, "Output path.", "path"});

    parser.addOption({{"fs", "fftsize"},
                      "Window Size Index (0 = 8192, 1 = 16384, 2 = 32768, 3 = 65536, 4 = 131072).",
                      "index",
                      "2"});

    parser.addOption({{"res", "resolution"}, "Resolution/Padding (1, 2, 4, 8).", "factor", "1"});
    parser.addOption(
        {{"to", "timeoverlap"}, "Time Overlap Divisor (4, 8, 16, 32, 64, 128).", "factor", "4"});

    parser.addOption({{"t", "threshold"}, "Spectral Gate Threshold (dB).", "db", "-60"});

    parser.addOption({"ampMin", "Make signal type contrast (dB).", "db", "-90"});
    parser.addOption({"ampMax", "Make signal type contrast (dB).", "db", "-60"});

    parser.addOption({"format", "wav, flac, aiff.", "fmt", "wav"});
    parser.addOption({"bit", "16, 24, 32.", "depth", "24"});

    parser.process(app);

    ProcessSettings settings;
    settings.inputPath = parser.value("i");
    settings.outputPath = parser.value("o");

    if (settings.inputPath.isEmpty() || settings.outputPath.isEmpty())
        return 1;

    QDir outDir(settings.outputPath);
    if (!outDir.exists())
        outDir.mkpath(".");

    int fsIdx = parser.value("fs").toInt();
    switch (fsIdx) {
    case 0:
        settings.windowSize = 8192;
        break;
    case 1:
        settings.windowSize = 16384;
        break;
    case 2:
        settings.windowSize = 32768;
        break;
    case 3:
        settings.windowSize = 65536;
        break;
    case 4:
        settings.windowSize = 131072;
        break;
    default:
        settings.windowSize = 32768;
        break;
    }

    settings.resolutionFactor = parser.value("res").toInt();
    if (settings.resolutionFactor != 1 && settings.resolutionFactor != 2
        && settings.resolutionFactor != 4 && settings.resolutionFactor != 8) {
        if (settings.resolutionFactor > 8)
            settings.resolutionFactor = 8;
        else if (settings.resolutionFactor < 1)
            settings.resolutionFactor = 1;
    }

    settings.hopDivisor = parser.value("to").toInt();
    if (settings.hopDivisor < 2)
        settings.hopDivisor = 4;

    settings.thresholdDB = parser.value("t").toDouble();

    QString fmtStr = parser.value("format").toLower();
    if (fmtStr == "flac")
        settings.outputFormat = SF_FORMAT_FLAC;
    else if (fmtStr == "aiff")
        settings.outputFormat = SF_FORMAT_AIFF;
    else
        settings.outputFormat = SF_FORMAT_WAV;

    QString bitStr = parser.value("bit");
    if (bitStr == "24")
        settings.outputBitDepth = SF_FORMAT_PCM_24;
    else if (bitStr == "32")
        settings.outputBitDepth = SF_FORMAT_PCM_32;
    else if (bitStr == "float")
        settings.outputBitDepth = SF_FORMAT_FLOAT;
    else
        settings.outputBitDepth = SF_FORMAT_PCM_16;

    QFileInfo inInfo(settings.inputPath);
    QList<FileJob> jobs;

    if (inInfo.isFile()) {
        FileJob job;
        job.inputPath = settings.inputPath;
        job.settings = settings;
        jobs.append(job);
    } else if (inInfo.isDir()) {
        QDirIterator it(settings.inputPath, {"*.wav", "*.aiff", "*.flac"}, QDir::Files);
        while (it.hasNext()) {
            FileJob job;
            job.inputPath = it.next();
            job.settings = settings;
            jobs.append(job);
        }
    }

    QtConcurrent::blockingMap(jobs, processFileJob);

    return 0;
}
