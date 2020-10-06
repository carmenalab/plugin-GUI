/*
------------------------------------------------------------------

This file is part of the Open Ephys GUI
Copyright (C) 2013 Open Ephys

------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __SPIKESORTBOXES_H
#define __SPIKESORTBOXES_H

#include "SpikeSorterEditor.h"
#include <algorithm>    // std::sort
#include <list>
#include <queue>
#include <atomic>
#include <utility>
#include <sstream>
#include <nlohmann/json.hpp>
#include <unordered_set>

using json = nlohmann::json;

class SorterSpikeContainer : public ReferenceCountedObject
{
public:
	//This invalidates the original SpikeEventPtr, so be careful
	SorterSpikeContainer(const SpikeChannel* channel, SpikeEvent::SpikeBuffer& data, int64 timestamp);
	SorterSpikeContainer() = delete;

	const float* getData() const;
	const SpikeChannel* getChannel() const;
	int64 getTimestamp() const;
	uint8 color[3];
	std::vector<float> pcProj;
	uint16 sortedId;
private:
	int64 timestamp;
	HeapBlock<float> data;
	const SpikeChannel* chan;
};
typedef ReferenceCountedObjectPtr<SorterSpikeContainer> SorterSpikePtr;
typedef ReferenceCountedArray<SorterSpikeContainer, CriticalSection> SorterSpikeArray;

class PCAcomputingThread;
class UniqueIDgenerator;
class PointD
{
public:

    PointD();
    PointD(float x, float y);
    PointD(const std::vector<float>& dims);
    PointD(const PointD& P);
    const PointD operator+(const PointD& c1) const;
    PointD& operator+=(const PointD& rhs);
    PointD& operator-=(const PointD& rhs);

    const PointD operator-(const PointD& c1) const;
    const PointD operator*(const PointD& c1) const;

    float cross(PointD c) const;

    float operator[](int index) const;

    void set(const std::initializer_list<double>& dims) {
        dims_.clear();
        for (const double &dim : dims) {
            dims_.push_back((float) dim);
        }
    }

    void set(const std::vector<float>& dims) {
        dims_ = dims;
    }

    std::vector<float> dims() const {
        return dims_;
    }
private:
    std::vector<float> dims_;
};


class Box
{
public:
    Box();
    Box(int channel);
    Box(float X, float Y, float W, float H, int ch=0);
    bool LineSegmentIntersection(PointD p11, PointD p12, PointD p21, PointD p22);
    bool isWaveFormInside(SorterSpikePtr so);
    double x,y,w,h; // x&w and specified in microseconds. y&h in microvolts
    int channel;
};


/************************/
class Histogram
{
public:
    Histogram();
    Histogram(int N, double T0, double T1);
    ~Histogram();

    void setParameters(int N, double T0, double T1);
    std::vector<int> getCounter();
    void reset();
    void update(double x);

    int Max;
    double t0, t1;
    std::vector<double> Time;
    int numBins;
    std::vector<int> Counter;

};

class RunningStats
{
public:
    RunningStats();
    ~RunningStats();
    void resizeWaveform(int newlength);
    void reset();
    Histogram getHistogram();
    std::vector<double> getMean(int index);
    std::vector<double> getStandardDeviation(int index);
    void update(SorterSpikePtr so);
    bool queryNewData();

    double LastSpikeTime;
    bool newData;
    Histogram hist;
    std::vector<std::vector<double> > WaveFormMean,WaveFormSk,WaveFormMk;
    double numSamples;


};

// Box unit defines a single unit (with multiple boxes)
// Each box can be on a different channel
class BoxUnit
{
public:
    BoxUnit();
    BoxUnit(int ID, int localID);
    BoxUnit(Box B, int ID, int localID);
    bool isWaveFormInsideAllBoxes(SorterSpikePtr so);
    bool isActivated();
    void activateUnit();
    void deactivateUnit();
    double getNumSecondsActive();
    void toggleActive();
    void addBox(Box b);
    void addBox();
    int getNumBoxes();
    void modifyBox(int boxindex, Box b);
    bool deleteBox(int boxindex);
    Box getBox(int box);
    void setBox(int boxid, Box B);
    void setBoxPos(int boxid, PointD P);
    void setBoxSize(int boxid, double W, double H);
    void MoveBox(int boxid, int dx, int dy);
    std::vector<Box> getBoxes();
    int getUnitID();
    int getLocalID();
	void updateWaveform(SorterSpikePtr so);
    static void setDefaultColors(uint8_t col[3], int ID);
    void resizeWaveform(int newlength);
public:
    int UnitID;
    int localID; // used internally, for colors and position.
    std::vector<Box> lstBoxes;
    uint8_t ColorRGB[3];
    RunningStats WaveformStat;
    bool Active;
    juce::int64 Activated_TS_S;
    Time timer;

};

/*
class PCAjob
{
public:
PCAjob();
};*/


class PCAResults {
public:

    PCAResults() : PCAResults(0) {}

    PCAResults(int num_pcs) : num_pcs_(num_pcs) {}

    PCAResults(
            int num_pcs,
            std::vector<std::vector<float>>  pcs,
            std::vector<float>  mean,
            std::vector<float>  pc_mins,
            std::vector<float>  pc_maxes)
            : num_pcs_(num_pcs),
              pcs_(std::move(pcs)),
              mean_(std::move(mean)),
              pc_mins_(std::move(pc_mins)),
              pc_maxes_(std::move(pc_maxes)) {}

    int num_pcs() const {
        return num_pcs_;
    }

    const std::vector<std::vector<float>> &pcs() const {
        return pcs_;
    }

    const std::vector<float> &mean() const {
        return mean_;
    }

    const std::vector<float> &pc_mins() const {
        return pc_mins_;
    }

    const std::vector<float> &pc_maxes() const {
        return pc_maxes_;
    }

    bool is_populated() const {
        return num_pcs_ > 0 &&
               pcs_.size() == num_pcs_ &&
               pc_mins_.size() == num_pcs_ &&
               pc_maxes_.size() == num_pcs_;
    }

    void clear() {
        pcs_.clear();
        pc_mins_.clear();
        pc_maxes_.clear();
        mean_.clear();
    }

    json ToJson() const;
    bool FromJson(const json& value, int expected_pc_vector_samples);

    explicit operator std::string() const {
        std::stringstream ss;
        ss << "PCAResults[num_pcs=" << num_pcs_ << ",is_populated=" << is_populated() << "]";
        return ss.str();
    }

private:
    int num_pcs_;
    std::vector<std::vector<float>> pcs_;
    std::vector<float> mean_;
    std::vector<float> pc_mins_;
    std::vector<float> pc_maxes_;
};

class PCAjob : public ReferenceCountedObject {
public:
    PCAjob(SorterSpikeArray &_spikes,
           PCAResults *results,
           std::atomic<bool> &_reportDone);

    ~PCAjob();

    void computeCov();

    void computeSVD();

    float** cov;
    SorterSpikeArray spikes;
    std::atomic<bool>& reportDone;

private:
    int svdcmp(float** a, int nRows, int nCols, float* w, float** v);
    float pythag(float a, float b);
    int dim;
    PCAResults *results_;
};

typedef ReferenceCountedObjectPtr<PCAjob> PCAJobPtr;
typedef ReferenceCountedArray<PCAjob, CriticalSection> PCAJobArray;

class cShape
{
public:
    virtual ~cShape() = default;
    virtual bool isPointInside(PointD p) const = 0;
    virtual std::vector<PointD>& pts() = 0;
    virtual PointD& offset() = 0;
    virtual bool CanSerializeToXml() const = 0;
    virtual json ToJson() const = 0;
    virtual bool FromJson(const json& value) = 0;
};


class cPolygon : public cShape
{
public:
    bool isPointInside(PointD p) const override;
    json ToJson() const override;
    bool FromJson(const json& value) override;

    std::vector<PointD>& pts() override {
        return pts_;
    }
    PointD& offset() override {
        return offset_;
    };

    bool CanSerializeToXml() const override {
        return true;
    }
private:
    std::vector<PointD> pts_;
    PointD offset_;
};

class cEllipse : public cShape
{
public:
    explicit cEllipse() {
        n_dims = 0;
    };

    cEllipse(const PointD& center,
             std::vector<std::vector<float>>  rotation,
             const std::vector<float>& radii);

    bool isPointInside(PointD p) const override;
    json ToJson() const override;
    bool FromJson(const json& value) override;

    std::vector<PointD>& pts() override {
        return sample_pts_;
    }

    PointD& offset() override {
        return center_;
    };

    bool CanSerializeToXml() const override {
        // Can be done, just not needed at the moment (always loaded via the API for now)
        return false;
    }
private:
    int n_dims;
    PointD center_;
    std::vector<std::vector<float>> rotation_;
    std::vector<float> radii_;
    std::vector<PointD> sample_pts_;

    void populate_sample_pts();
};



class PCAcomputingThread : juce::Thread
{
public:
    PCAcomputingThread();
    void run(); // computes PCA on waveforms
    void addPCAjob(PCAJobPtr job);

private:
    PCAJobArray jobs;
	CriticalSection lock;
};

class PCAUnit
{
public:
    PCAUnit();
    PCAUnit(int ID, int localID);
    ~PCAUnit();
    int getUnitID();
    int getLocalID();
	bool isWaveFormInsidePolygon(SorterSpikePtr so);
    bool isPointInsidePolygon(PointD p);
    void resizeWaveform(int newlength);
	void updateWaveform(SorterSpikePtr so);

    json ToJson() const;
    bool FromJson(const json& value);
public:
    int UnitID;
    int localID; // used internally, for colors and position.
    std::shared_ptr<cShape> poly;
    uint8_t ColorRGB[3];
    RunningStats WaveformStat;
    bool Active;
    juce::int64 Activated_TS_S;
    Time timer;
};

class PCASorting {
public:
    PCAResults results;
    std::vector<PCAUnit> units;

    Value ToValue() const;
    bool FromValue(const Value& value, int expected_pc_vector_samples);
};

// Sort spikes from a single electrode (which could have any number of channels)
// using the box method. Any electrode could have an arbitrary number of units specified.
// Each unit is defined by a set of boxes, which can be placed on any of the given channels.
class SpikeSortBoxes : Parameter::Listener
{
public:
    SpikeSortBoxes(UniqueIDgenerator *uniqueIDgenerator_,
                   PCAcomputingThread *pth,
                   int numch,
                   double SamplingRate,
                   int WaveFormLength,
                   Parameter *pca_parameter);
    ~SpikeSortBoxes();

    void resizeWaveform(int numSamples);

    // Parameter::Listener
    void parameterValueChanged(Value &valueThatWasChanged, const String &parameterName) override;

    // Populate the SorterSpikePtr's pcProj array if there are possible PC units.
    void projectOnPrincipalComponents(SorterSpikePtr so);

    // Try to sort the spike, populating the SorterSpikePtr's sortedID and color to the appropriate unit
    // if sorted.
	bool sortSpike(SorterSpikePtr so, bool PCAfirst);
    void RePCA();

    void DisablePeriodicPCA();
    void EnablePeriodicPCA();

    void addPCAunit(PCAUnit unit);
    int addBoxUnit(int channel);
    int addBoxUnit(int channel, Box B);

    void getPCArange(int pc_idx, float &min, float &max);
    void setPCArange(int pc_idx, float min, float max);

    void resetJobStatus();
    bool isPCAfinished();

    bool removeUnit(int unitID);

    void removeAllUnits();
    bool addBoxToUnit(int channel, int unitID);
    bool addBoxToUnit(int channel, int unitID, Box B);
    bool removeBoxFromUnit(int unitID, int boxIndex);
    int getNumBoxes(int unitID);
    std::vector<Box> getUnitBoxes(int unitID);
    std::vector<BoxUnit> getBoxUnits();
    std::vector<PCAUnit> getPCAUnits();

    void getUnitColor(int UnitID, uint8& R, uint8& G, uint8& B);
    void updateBoxUnits(std::vector<BoxUnit> _units);
    void updatePCAUnits(std::vector<PCAUnit> _units);
    int generateUnitID();
    int generateLocalID();
    void generateNewIDs();
    void setSelectedUnitAndBox(int unitID, int boxID);
    void getSelectedUnitAndBox(int& unitID, int& boxid);
    void saveCustomParametersToXml(XmlElement* electrodeNode);
    void loadCustomParametersFromXml(XmlElement* electrodeNode);
private:
    void synchronizePCAParameter();
    UniqueIDgenerator* uniqueIDgenerator;
    int numChannels, waveformLength;
    int selectedUnit, selectedBox;
    CriticalSection mut;
    std::vector<BoxUnit> boxUnits;

    PCASorting pca_sorting_;
    // Updated by the PCAcomputingThread; the process thread will deal with replacing the
    // "active" pca sorting instance.
    PCASorting replacement_pca_sorting_;

    SorterSpikeArray spikeBuffer;
    int bufferSize,spikeBufferIndex;
    PCAcomputingThread* computingThread;
    bool bPCAJobSubmitted,bRePCA, bShouldPCA;
    std::atomic<bool> bPCAjobFinished ;

    Parameter *pca_parameter_;
};

//Those are legacy methods from the old spike system that are must likely not needed in the new one
float spikeDataBinToMicrovolts(SorterSpikePtr  s, int bin, int ch = 0);
float spikeDataIndexToMicrovolts(SorterSpikePtr s, int index = 0);
float spikeTimeBinToMicrosecond(SorterSpikePtr s, int bin, int ch = 0);
int microSecondsToSpikeTimeBin(SorterSpikePtr s, float t, int ch = 0);


#endif // __SPIKESORTBOXES_H
