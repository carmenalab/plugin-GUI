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

#include <stdio.h>
#include <algorithm>
#include "SpikeSortBoxes.h"
#include "SpikeSorter.h"
#include <nlohmann/json.hpp>
#include <cmath>

using json = nlohmann::json;

PointD::PointD() : PointD(0.0, 0.0) {}
PointD::PointD(float x, float y) : PointD(std::vector<float>({x, y})) {}
PointD::PointD(const PointD& P) : PointD(P.dims_) {}
PointD::PointD(const std::vector<float> &dims) {
    dims_ = dims;
}

PointD& PointD::operator+=(const PointD& rhs)
{
    jassert(dims_.size() == rhs.dims_.size());
    for (int i = 0; i < dims_.size(); i++) {
        dims_[i] += rhs.dims_[i];
    }
    return *this;
}

PointD& PointD::operator-=(const PointD& rhs)
{
    for (int i = 0; i < dims_.size(); i++) {
        dims_[i] -= rhs.dims_[i];
    }
    return *this;
}

const PointD PointD::operator+(const PointD& other) const
{
    PointD result = *this;
    result += other;
    return result;
}


const PointD PointD::operator-(const PointD& other) const
{
    PointD result = *this;
    result -= other;
    return result;
}


const PointD PointD::operator*(const PointD& other) const
{
    PointD result = *this;
    for (int i = 0; i < result.dims_.size(); i++) {
        result.dims_[i] *= other.dims_[i];
    }
    return result;

}

float PointD::cross(PointD c) const
{
    jassert(dims_.size() == c.dims_.size());
    jassert(dims_.size() == 2);
    return dims_[0] * c[1] - dims_[1] * c[0];
}

float PointD::operator[](int index) const {
    return dims_[index];
}

/**************************************/


Box::Box()
{
    x = -0.2; // in ms
    y = -10; // in uV
    w = 0.5; // in ms
    h = 70; // in uV
    channel=0;
}


Box::Box(int ch)
{
    x = -0.2; // in ms
    y = -10; // in uV
    w = 0.5; // in ms
    h = 70; // in uV
    channel = ch;
}

Box::Box(float X, float Y, float W, float H, int ch)
{
    x = X;
    y = Y;
    w = W;
    h = H;
    channel = ch;
}

bool Box::LineSegmentIntersection(PointD p11, PointD p12, PointD p21, PointD p22)
{
    PointD r = (p12 - p11);
    PointD s = (p22 - p21);
    PointD q = p21;
    PointD p = p11;
    double rs = r.cross(s);
    double eps = 1e-6;
    if (fabs(rs) < eps)
        return false; // lines are parallel
    double t = (q - p).cross(s) / rs;
    double u = (q - p).cross(r) / rs;
    return (t>=0&&t<=1 &&u>0&&u<=1);
}

#ifndef MAX
#define MAX(x,y)((x)>(y))?(x):(y)
#endif

#ifndef MIN
#define MIN(x,y)((x)<(y))?(x):(y)
#endif

bool Box::isWaveFormInside(SorterSpikePtr so)
{
    PointD BoxTopLeft(x, y);
    PointD BoxBottomLeft(x, (y - h));

    PointD BoxTopRight(x + w, y);
    PointD BoxBottomRight(x + w, (y - h));

    // y,and h are given in micro volts.
    // x and w and given in micro seconds.

    // no point testing all wave form points. Just ones that are between x and x+w...
    int BinLeft = microSecondsToSpikeTimeBin(so,x);
    int BinRight = microSecondsToSpikeTimeBin(so,x+w);

    /*
    float minValue=1e10, maxValue=1e-10;
    for (int pt = 0; pt < so->nSamples; pt++)
    {
    	float v = spikeDataBinToMicrovolts(so, pt, channel);
    	minValue = MIN(minValue,v);
    	maxValue = MAX(maxValue,v);
    }
    */

    for (int pt = BinLeft; pt < BinRight; pt++)
    {
        PointD Pwave1(spikeTimeBinToMicrosecond(so,pt),spikeDataBinToMicrovolts(so, pt, channel));
        PointD Pwave2(spikeTimeBinToMicrosecond(so,pt+1),spikeDataBinToMicrovolts(so, pt+1, channel));

        bool bLeft = LineSegmentIntersection(Pwave1,Pwave2,BoxTopLeft,BoxBottomLeft) ;
        bool bRight = LineSegmentIntersection(Pwave1,Pwave2,BoxTopRight,BoxBottomRight);
        bool bTop = LineSegmentIntersection(Pwave1,Pwave2,BoxTopLeft,BoxTopRight);
        bool bBottom = LineSegmentIntersection(Pwave1, Pwave2, BoxBottomLeft, BoxBottomRight);
        if (bLeft || bRight || bTop || bBottom)
        {
            return true;
        }

    }
    return false;
}


BoxUnit::BoxUnit()
{

}
/**************************************/
void BoxUnit::setDefaultColors(uint8_t col[3], int ID)
{
    int IDmodule = (ID-1) % 6; // ID can't be zero
    const int colors[6][3] =
    {
        {0xFF,0xFF,0x00},
        {0x00,0xFF,0x00},
        {0x00, 0xFF, 0xFF},
        {0xFF, 0x00, 0x00},
        {0x00,0x00,0xFF},
        {0xFF,0x00,0xFF}
    };
    col[0] = colors[IDmodule][0];
    col[1] = colors[IDmodule][1];
    col[2] = colors[IDmodule][2];
}

BoxUnit::BoxUnit(int ID, int localID_): UnitID(ID), localID(localID_)
{
    std::cout << "Adding new box unit." << std::endl;
    Active = false;
    Activated_TS_S = -1;
    setDefaultColors(ColorRGB, localID);
    Box B(50, -20 - localID*20, 300, 40);
    addBox(B);
}

void BoxUnit::resizeWaveform(int newlength)
{
    WaveformStat.resizeWaveform(newlength);
}

BoxUnit::BoxUnit(Box B, int ID, int localID_) : UnitID(ID), localID(localID_)
{
    addBox(B);
}

bool BoxUnit::isWaveFormInsideAllBoxes(SorterSpikePtr so)
{
    for (int k=0; k< lstBoxes.size(); k++)
    {
        if (!lstBoxes[k].isWaveFormInside(so))
            return false;
    }
    return lstBoxes.size() == 0 ? false : true;
}

bool BoxUnit::isActivated()
{
    return Active;
}

void BoxUnit::activateUnit()
{
    Active = true;
    Activated_TS_S = timer.getHighResolutionTicks();
}

void BoxUnit::deactivateUnit()
{
    Active = false;
    Activated_TS_S = timer.getHighResolutionTicks();

}

double BoxUnit::getNumSecondsActive()
{
    if (!Active)
        return 0;
    else
        return (timer.getHighResolutionTicks() - Activated_TS_S) / timer.getHighResolutionTicksPerSecond();
}

void BoxUnit::toggleActive()
{
    if (Active)
        deactivateUnit();
    else
        activateUnit();
}

void BoxUnit::addBox(Box b)
{
    lstBoxes.push_back(b);
}

void BoxUnit::addBox()
{
    Box B(50 + 350 * lstBoxes.size(), -20 - UnitID * 20, 300, 40);
    lstBoxes.push_back(B);
}

int BoxUnit::getNumBoxes()
{
    return (int) lstBoxes.size();
}

void BoxUnit::modifyBox(int boxindex, Box b)
{
    lstBoxes[boxindex] = b;
}


bool BoxUnit::deleteBox(int boxindex)
{

    if (lstBoxes.size() > boxindex)
    {
        lstBoxes.erase(lstBoxes.begin()+boxindex);
        return true;
    }
    return false;
}

Box BoxUnit::getBox(int box)
{
    return lstBoxes[box];
}

void BoxUnit::setBox(int boxid, Box B)
{
    lstBoxes[boxid].x = B.x;
    lstBoxes[boxid].y = B.y;
    lstBoxes[boxid].w = B.w;
    lstBoxes[boxid].h = B.h;
}


void BoxUnit::setBoxPos(int boxid, PointD P)
{
    lstBoxes[boxid].x = P[0];
    lstBoxes[boxid].y = P[1];
}

void BoxUnit::setBoxSize(int boxid, double W, double H)
{
    lstBoxes[boxid].w = W;
    lstBoxes[boxid].h = H;
}

void BoxUnit::MoveBox(int boxid, int dx, int dy)
{
    lstBoxes[boxid].x += dx;
    lstBoxes[boxid].y += dy;
}

std::vector<Box> BoxUnit::getBoxes()
{
    return lstBoxes;
}

// Members
int BoxUnit::getUnitID()
{
    return UnitID;
}

int BoxUnit::getLocalID()
{
    return localID;
}


/************************/

Histogram::~Histogram()
{
    Time.clear();
    Counter.clear();
}

Histogram::Histogram()
{

}

void Histogram::setParameters(int N, double T0, double T1)
{

    t0 = T0;
    t1 = T1;
    numBins = N;
    Time.resize(N);
    Counter.resize(N);
    for (int k = 0; k < N; k++)
    {
        Time[k] = (double)k / (N - 1) * (T1 - T0) + T0;
        Counter[k] = 0;
    }
}

Histogram::Histogram(int N, double T0, double T1)
{
    t0 = T0;
    t1 = T1;
    numBins = N;
    Time.resize(N);
    Counter.resize(N);
    for (int k = 0; k < N; k++)
    {
        Time[k] = (double)k / (N - 1) * (T1 - T0) + T0;
        Counter[k] = 0;
    }
}

void Histogram::update(double x)
{
    int Bin = ((x - t0) / (t1 - t0) * (numBins-1));
    if (Bin >= 0 && Bin < numBins)
    {
        Counter[Bin]++;
        if (Counter[Bin] > Max)
        {
            Max = Counter[Bin];
        }
    }
}

void Histogram::reset()
{
    Max = 0;
    for (int k = 0; k < numBins; k++)
        Counter[k] = 0;
}

std::vector<int> Histogram::getCounter()
{
    return Counter;
}


/**********************/
// computes statistics about the unit (non-trial related statistics)



// Running variance...
//Mk = Mk-1+ (xk - Mk-1)/k
//Sk = Sk-1 + (xk - Mk-1)*(xk - Mk).
//For 2 ≤ k ≤ n, the kth estimate of the variance is s2 = Sk/(k - 1).
RunningStats::~RunningStats()
{
}
RunningStats::RunningStats()
{
    hist.setParameters(101, 0, 100); // Inter spike histogram. Fixed range [0..100 ms]
    numSamples = 0;
}

void RunningStats::reset()
{
    numSamples = 0;
    hist.reset();
}

Histogram RunningStats::getHistogram()
{
    return hist;
}

std::vector<double> RunningStats::getMean(int index)
{
    std::vector<double> m;

    if (numSamples == 0)
    {
        return m;
    }

    int numSamplesInWaveForm = (int) WaveFormMean[0].size();
    m.resize(numSamplesInWaveForm);

    for (int k = 0; k < numSamplesInWaveForm; k++)
        m[k] = WaveFormMean[index][k];
    return m;
}

std::vector<double> RunningStats::getStandardDeviation(int index)
{
    std::vector<double> WaveFormVar;

    if (numSamples == 0)
    {
        return WaveFormVar;
    }
    int numSamplesInWaveForm = (int) WaveFormMean[0].size();
    WaveFormVar.resize(numSamplesInWaveForm);

    for (int j = 0; j < numSamplesInWaveForm; j++)
    {
        if (numSamples - 1 == 0)
            WaveFormVar[j] = 0;
        else
            WaveFormVar[j] = sqrt(WaveFormSk[index][j] / (numSamples - 1));
    }
    return WaveFormVar;
}


void RunningStats::resizeWaveform(int newlength)
{
    numSamples = 0; // this should ensure that update reallocates upon the next update.
}

void RunningStats::update(SorterSpikePtr so)
{
    double ts = so->getTimestamp()/so->getChannel()->getSampleRate();
    if (numSamples == 0)
    {
        LastSpikeTime = ts;
    }
    else
    {
        hist.update(1000.0 * (ts - LastSpikeTime));
        LastSpikeTime = ts;
    }

    newData = true;
	int nChannels = so->getChannel()->getNumChannels();
	int nSamples = so->getChannel()->getTotalSamples();
    if (numSamples == 0)
    {
        // allocate
        WaveFormMean.resize(nChannels);
        WaveFormSk.resize(nChannels);
        WaveFormMk.resize(nChannels);
        for (int k=0; k<nChannels; k++)
        {
            WaveFormMean[k].resize(nSamples);
            WaveFormSk[k].resize(nSamples);
            WaveFormMk[k].resize(nSamples);
        }

        for (int i = 0; i < nChannels; i++)
        {
            for (int j = 0; j < nSamples; j++)
            {
				WaveFormMean[i][j] = so->getData()[j + i*nSamples];
                WaveFormSk[i][j] = 0;
				WaveFormMk[i][j] = so->getData()[j + i*nSamples];
            }
        }
        numSamples += 1.0F;
        return;
    }
    // running mean
    for (int i = 0; i < nChannels; i++)
    {
        for (int j = 0; j < nSamples; j++)
        {
			WaveFormMean[i][j] = (numSamples * WaveFormMean[i][j] + so->getData()[j + i*nSamples]) / (numSamples + 1);
			WaveFormMk[i][j] += (so->getData()[j + i*nSamples] - WaveFormMk[i][j]) / numSamples;
			WaveFormSk[i][j] += (so->getData()[j + i*nSamples] - WaveFormMk[i][j]) * (so->getData()[j + i*nSamples] - WaveFormMk[i][j]);
        }
    }
    numSamples += 1.0F;
}


bool RunningStats::queryNewData()
{
    if (newData == false)
        return false;
    newData = false;
    return true;
}

void BoxUnit::updateWaveform(SorterSpikePtr so)
{
    WaveformStat.update(so);
}


/*
        public bool QueryNewData()
        {
            return WaveFormStat.QueryNewData();
        }

        public Histogram GetInterSpikeHistogram()
        {
            return WaveFormStat.GetHistogram();
        }

        public void GetWaveFormMeanStd(out double[] Mean, out double[] Std, int index)
        {
            Mean = WaveFormStat.GetMean(index);
            Std = WaveFormStat.GetStandardDeviation(index);
        }

        public void ResetWaveForm()
        {
            WaveFormStat.Reset();
        }


    }
	*/

/***********************************************/

SpikeSortBoxes::SpikeSortBoxes(UniqueIDgenerator *uniqueIDgenerator_,
                               PCAcomputingThread *pth,
                               int numch,
                               double SamplingRate,
                               int WaveFormLength,
                               Parameter *pca_parameter) {
    uniqueIDgenerator = uniqueIDgenerator_;
    computingThread = pth;
    bufferSize = 200;
    spikeBufferIndex = -1;
    bPCAcomputed = false;
    bPCAJobSubmitted = false;
    bPCAjobFinished = false;
    selectedUnit = -1;
    selectedBox = -1;
    bRePCA = false;
    numChannels = numch;
    waveformLength = WaveFormLength;

    pca_parameter_ = pca_parameter;

    // Sane default for now.
    pca_sorting_.results.num_pcs = 2;

    for (int n = 0; n < bufferSize; n++)
    {
        spikeBuffer.add(nullptr);
    }

    pca_parameter_->addListener(this);
    synchronizePCAParameter();
}

void SpikeSortBoxes::synchronizePCAParameter() {
    pca_parameter_->setValue(pca_sorting_.ToValue(), 0);
}

void SpikeSortBoxes::parameterValueChanged(Value &valueThatWasChanged) {
    if (pca_parameter_->getNumChannels() != 1) {
        // Something invalid? Only one channel supported.
        return;
    }

    juce::String parameter_value_str = pca_parameter_->getValue(0).toString();

    PCASorting new_sorting;
    bool success = new_sorting.FromValue(Value(pca_parameter_->getValue(0)), numChannels * waveformLength);
    if (success) {
        pca_sorting_ = new_sorting;
    }
    // In case the parameter was invalid or was otherwise unable to construct sorting, ensure the parameter is in
    // sync with the "real".
    synchronizePCAParameter();
}

void SpikeSortBoxes::resizeWaveform(int numSamples)
{
    const ScopedLock myScopedLock(mut);

    waveformLength = numSamples;
    pca_sorting_.results.clear();
    spikeBuffer.clear();
    for (int n = 0; n < bufferSize; n++)
    {
        spikeBuffer.add(nullptr);
    }
    bPCAcomputed = false;
    spikeBufferIndex = -1;
	bPCAJobSubmitted = false;
	bPCAjobFinished = false;
	selectedUnit = -1;
	selectedBox = -1;
	bRePCA = false;

	for (auto &pcaUnit : pca_sorting_.units) {
	    pcaUnit.resizeWaveform(waveformLength);
	}
    for (int k=0; k<boxUnits.size(); k++)
    {
        boxUnits[k].resizeWaveform(waveformLength);
    }

    synchronizePCAParameter();
}



void SpikeSortBoxes::loadCustomParametersFromXml(XmlElement* electrodeNode)
{

    forEachXmlChildElement(*electrodeNode, spikesortNode)
    {
        if (spikesortNode->hasTagName("SPIKESORTING"))
        {
            //int numBoxUnit  = spikesortNode->getIntAttribute("numBoxUnits");
            //int numPCAUnit  = spikesortNode->getIntAttribute("numPCAUnits");
            selectedUnit  = spikesortNode->getIntAttribute("selectedUnit");
            selectedBox =  spikesortNode->getIntAttribute("selectedBox");

            pca_sorting_ = PCASorting();
            PCAResults &pca_results = pca_sorting_.results;
            std::vector<PCAUnit> &pcaUnits = pca_sorting_.units;

            boxUnits.clear();

            forEachXmlChildElement(*spikesortNode, UnitNode)
            {
                if (UnitNode->hasTagName("PCA"))
                {
                    numChannels = UnitNode->getIntAttribute("numChannels");
                    waveformLength = UnitNode->getIntAttribute("waveformLength");
                    bPCAjobFinished = UnitNode->getBoolAttribute("PCAjobFinished");
                    bPCAcomputed = UnitNode->getBoolAttribute("PCAcomputed");

                    pca_results.num_pcs = UnitNode->getIntAttribute("numPCs", 2);

                    bool is_old_style_pcs = false;
                    if (UnitNode->hasAttribute("pc1min")) {
                        // "Old" style with 2 hardcoded PCs, detected by presence of "old" tag
                        jassert(pca_results.num_pcs == 2);
                        is_old_style_pcs = true;
                        pca_results.pc_mins.push_back(UnitNode->getDoubleAttribute("pc1min"));
                        pca_results.pc_mins.push_back(UnitNode->getDoubleAttribute("pc2min"));

                        pca_results.pc_maxes.push_back(UnitNode->getDoubleAttribute("pc1max"));
                        pca_results.pc_maxes.push_back(UnitNode->getDoubleAttribute("pc2max"));
                    }

                    forEachXmlChildElement(*UnitNode, dimNode)
                    {

                        if (dimNode->hasTagName("PCA_DIM"))
                        {
                            // "Old" style with 2 hardcoded PCs if it has "PCA_DIM"
                            jassert(is_old_style_pcs);
                            if (pca_results.pcs.empty()) {
                                pca_results.pcs.emplace_back();
                                pca_results.pcs.emplace_back();
                            }
                            pca_results.pcs[0].push_back(dimNode->getDoubleAttribute("pc1"));
                            pca_results.pcs[1].push_back(dimNode->getDoubleAttribute("pc2"));
                        } else if (dimNode->hasTagName("PC")) {
                            jassert(!is_old_style_pcs);
                            pca_results.pc_maxes.push_back(dimNode->getDoubleAttribute("max_projection"));
                            pca_results.pc_mins.push_back(dimNode->getDoubleAttribute("min_projection"));
                            std::vector<float> pc_vector;
                            forEachXmlChildElement(*dimNode, pc_element_node)
                            {
                                if (pc_element_node->hasTagName("PC_ELEMENT")) {
                                    pc_vector.push_back(pc_element_node->getDoubleAttribute("value"));
                                }
                            }
                            pca_results.pcs.push_back(pc_vector);
                        }
                    }

                    // If any of the dimensions across the various PC fields mismatch, then start fresh.
                    if (!pca_results.is_populated()) {
                        pca_results.clear();
                    }
                }

                if (UnitNode->hasTagName("BOXUNIT"))
                {
                    BoxUnit boxUnit;
                    boxUnit.UnitID = UnitNode->getIntAttribute("UnitID");
                    boxUnit.ColorRGB[0] = UnitNode->getIntAttribute("ColorR");
                    boxUnit.ColorRGB[1] = UnitNode->getIntAttribute("ColorG");
                    boxUnit.ColorRGB[2] = UnitNode->getIntAttribute("ColorB");
                    int numBoxes = UnitNode->getIntAttribute("NumBoxes");
                    boxUnit.lstBoxes.resize(numBoxes);
                    int boxCounter = 0;
                    forEachXmlChildElement(*UnitNode, boxNode)
                    {
                        if (boxNode->hasTagName("BOX"))
                        {
                            Box box;
                            box.channel = boxNode->getIntAttribute("ch");
                            box.x = boxNode->getDoubleAttribute("x");
                            box.y = boxNode->getDoubleAttribute("y");
                            box.w = boxNode->getDoubleAttribute("w");
                            box.h = boxNode->getDoubleAttribute("h");
                            boxUnit.lstBoxes[boxCounter++] = box;
                        }
                    }
                    // add box unit
                    boxUnits.push_back(boxUnit);
                }
                if (UnitNode->hasTagName("PCAUNIT"))
                {
                    PCAUnit pcaUnit;

                    pcaUnit.UnitID = UnitNode->getIntAttribute("UnitID");
                    pcaUnit.ColorRGB[0] = UnitNode->getIntAttribute("ColorR");
                    pcaUnit.ColorRGB[1] = UnitNode->getIntAttribute("ColorG");
                    pcaUnit.ColorRGB[2] = UnitNode->getIntAttribute("ColorB");

                    int numPolygonPoints = UnitNode->getIntAttribute("PolygonNumPoints");
                    pcaUnit.poly = std::make_shared<cPolygon>();
                    std::vector<PointD> &pts = pcaUnit.poly->pts();
                    pts.resize(numPolygonPoints);

                    PointD &offset = pcaUnit.poly->offset();
                    offset.set({UnitNode->getDoubleAttribute("PolygonOffsetX"),
                                UnitNode->getDoubleAttribute("PolygonOffsetY")});
                    // read polygon
                    int pointCounter = 0;
                    forEachXmlChildElement(*UnitNode, polygonPoint)
                    {
                        if (polygonPoint->hasTagName("POLYGON_POINT"))
                        {
                            pts[pointCounter].set({polygonPoint->getDoubleAttribute("pointX"),
                                                   polygonPoint->getDoubleAttribute("pointY")});
                            pointCounter++;
                        }
                    }
                    // add polygon unit
                    pcaUnits.push_back(pcaUnit);
                }
            }
        }
    }

    synchronizePCAParameter();
}

void SpikeSortBoxes::saveCustomParametersToXml(XmlElement* electrodeNode)
{

    const std::vector<PCAUnit> &pcaUnits = pca_sorting_.units;
    const PCAResults &results = pca_sorting_.results;

    XmlElement* spikesortNode = electrodeNode->createNewChildElement("SPIKESORTING");
    spikesortNode->setAttribute("numBoxUnits", (int)boxUnits.size());
    spikesortNode->setAttribute("numPCAUnits", (int)pcaUnits.size());
    spikesortNode->setAttribute("selectedUnit",selectedUnit);
    spikesortNode->setAttribute("selectedBox",selectedBox);


    XmlElement* pcaNode = spikesortNode->createNewChildElement("PCA");
    pcaNode->setAttribute("numChannels",numChannels);
    pcaNode->setAttribute("waveformLength",waveformLength);
    pcaNode->setAttribute("PCAjobFinished", bPCAjobFinished);
    pcaNode->setAttribute("PCAcomputed", bPCAcomputed);
    pcaNode->setAttribute("numPCs", results.num_pcs);

    if (results.is_populated()) {
        for (int pc_idx = 0; pc_idx < results.num_pcs; pc_idx++) {
            XmlElement* pc = pcaNode->createNewChildElement("PC");
            pc->setAttribute("min_projection", results.pc_mins[pc_idx]);
            pc->setAttribute("max_projection", results.pc_maxes[pc_idx]);

            auto pc_vector = results.pcs[pc_idx];
            for (float k : pc_vector) {
                XmlElement *pc_element = pc->createNewChildElement("PC_ELEMENT");
                pc_element->setAttribute("value", k);
            }
        }
    }

    for (int boxUnitIter=0; boxUnitIter<boxUnits.size(); boxUnitIter++)
    {
        XmlElement* BoxUnitNode = spikesortNode->createNewChildElement("BOXUNIT");

        BoxUnitNode->setAttribute("UnitID",boxUnits[boxUnitIter].UnitID);
        BoxUnitNode->setAttribute("ColorR",boxUnits[boxUnitIter].ColorRGB[0]);
        BoxUnitNode->setAttribute("ColorG",boxUnits[boxUnitIter].ColorRGB[1]);
        BoxUnitNode->setAttribute("ColorB",boxUnits[boxUnitIter].ColorRGB[2]);
        BoxUnitNode->setAttribute("NumBoxes", (int)boxUnits[boxUnitIter].lstBoxes.size());
        for (int boxIter=0; boxIter<boxUnits[boxUnitIter].lstBoxes.size(); boxIter++)
        {
            XmlElement* BoxNode = BoxUnitNode->createNewChildElement("BOX");
            BoxNode->setAttribute("ch", (int)boxUnits[boxUnitIter].lstBoxes[boxIter].channel);
            BoxNode->setAttribute("x", (int)boxUnits[boxUnitIter].lstBoxes[boxIter].x);
            BoxNode->setAttribute("y", (int)boxUnits[boxUnitIter].lstBoxes[boxIter].y);
            BoxNode->setAttribute("w", (int)boxUnits[boxUnitIter].lstBoxes[boxIter].w);
            BoxNode->setAttribute("h", (int)boxUnits[boxUnitIter].lstBoxes[boxIter].h);
        }
    }

    for (int pcaUnitIter=0; pcaUnitIter<pcaUnits.size(); pcaUnitIter++)
    {
        const PCAUnit &pcaUnit = pcaUnits[pcaUnitIter];
        if (!pcaUnit.poly->CanSerializeToXml()) {
            continue;
        }

        XmlElement* PcaUnitNode = spikesortNode->createNewChildElement("PCAUNIT");

        PcaUnitNode->setAttribute("UnitID", pcaUnit.UnitID);
        PcaUnitNode->setAttribute("ColorR", pcaUnit.ColorRGB[0]);
        PcaUnitNode->setAttribute("ColorG", pcaUnit.ColorRGB[1]);
        PcaUnitNode->setAttribute("ColorB", pcaUnit.ColorRGB[2]);
        std::vector<PointD> &pts = pcaUnit.poly->pts();
        PointD &offset = pcaUnit.poly->offset();

        PcaUnitNode->setAttribute("PolygonNumPoints", (int) pts.size());
        PcaUnitNode->setAttribute("PolygonOffsetX", (int) offset[0]);
        PcaUnitNode->setAttribute("PolygonOffsetY",(int) offset[1]);

        for (int p=0; p < pts.size(); p++)
        {
            XmlElement* PolygonNode = PcaUnitNode->createNewChildElement("POLYGON_POINT");
            PolygonNode->setAttribute("pointX", pts[p][0]);
            PolygonNode->setAttribute("pointY", pts[p][1]);
        }
    }

}

SpikeSortBoxes::~SpikeSortBoxes()
{
}

void SpikeSortBoxes::setSelectedUnitAndBox(int unitID, int boxID)
{
    selectedUnit = unitID;
    selectedBox = boxID;
}

void SpikeSortBoxes::getSelectedUnitAndBox(int& unitID, int& boxid)
{
    unitID = selectedUnit;
    boxid = selectedBox;
}

void SpikeSortBoxes::projectOnPrincipalComponents(SorterSpikePtr so)
{
    spikeBufferIndex++;
    spikeBufferIndex %= bufferSize;
    spikeBuffer.set(spikeBufferIndex, so);
    if (bPCAjobFinished)
    {
        // PCA results have been updated by the PCA job, so ensure the parameter is updated too:
        synchronizePCAParameter();
        bPCAcomputed = true;
    }

    if (bPCAcomputed)
    {
        so->pcProj.clear();
        int num_pcs = pca_sorting_.results.num_pcs;
        so->pcProj.resize(num_pcs, 0.0);
        for (int pc_idx = 0; pc_idx < pca_sorting_.results.num_pcs; pc_idx++) {
            for (int k=0; k<so->getChannel()->getNumChannels()*so->getChannel()->getTotalSamples(); k++)
            {
                float v = spikeDataIndexToMicrovolts(so, k);
                so->pcProj[pc_idx] += pca_sorting_.results.pcs[pc_idx][k]* v;
            }
        }
    }
    else
    {
        // add a spike object to the buffer.
        // if we have enough spikes, start the PCA computation thread.
        bool spikeBufferFull = spikeBufferIndex == bufferSize - 1;
        if (bShouldPCA && ((spikeBufferFull && !bPCAcomputed && !bPCAJobSubmitted) || bRePCA)) {
            bPCAJobSubmitted = true;
            bPCAcomputed = false;
            bRePCA = false;
            // submit a new job to compute the spike buffer.
            PCAJobPtr job = new PCAjob(spikeBuffer, &pca_sorting_.results, bPCAjobFinished);
            computingThread->addPCAjob(job);
        }
    }
}

void SpikeSortBoxes::getPCArange(int pc_idx, float& min,float& max)
{
    PCAResults &results = pca_sorting_.results;
    if (results.pc_mins.size() <= pc_idx || results.pc_maxes.size() <= pc_idx) {
        min = -1.0;
        max = 1.0;
    } else {
        min = results.pc_mins[pc_idx];
        max = results.pc_maxes[pc_idx];
    }
}

void SpikeSortBoxes::setPCArange(int pc_idx, float min,float max)
{
    PCAResults &results = pca_sorting_.results;
    if (results.pc_mins.size() <= pc_idx || results.pc_maxes.size() <= pc_idx) {
        return;
    }
    results.pc_mins[pc_idx] = min;
    results.pc_maxes[pc_idx] = max;
}

void SpikeSortBoxes::resetJobStatus()
{
    bPCAjobFinished = false;
}

bool SpikeSortBoxes::isPCAfinished()
{
    return bPCAjobFinished;
}
void SpikeSortBoxes::RePCA()
{
    bPCAcomputed = false;
    bPCAJobSubmitted = false;
    bRePCA = true;
}

void SpikeSortBoxes::DisablePeriodicPCA() {
    bRePCA = false;
    bShouldPCA = false;
}

void SpikeSortBoxes::EnablePeriodicPCA() {
    bRePCA = true;
    bShouldPCA = true;
}

void SpikeSortBoxes::addPCAunit(PCAUnit unit)
{
    const ScopedLock myScopedLock(mut);
    //StartCriticalSection();
    pca_sorting_.units.push_back(unit);
    synchronizePCAParameter();
    //EndCriticalSection();
}

// Adds a new unit with a single box at some default location.
// returns the unit id
int SpikeSortBoxes::addBoxUnit(int channel)
{
    const ScopedLock myScopedLock(mut);
    //StartCriticalSection();
    int unusedID = uniqueIDgenerator->generateUniqueID(); //generateUnitID();
    BoxUnit unit(unusedID, generateLocalID());
    boxUnits.push_back(unit);
    setSelectedUnitAndBox(unusedID, 0);
    //EndCriticalSection();
    return unusedID;
}
/*
void  SpikeSortBoxes::StartCriticalSection()
{
	mut.enter();
}

void  SpikeSortBoxes::EndCriticalSection()
{
	mut.exit();
}
*/
int SpikeSortBoxes::addBoxUnit(int channel, Box B)
{
    const ScopedLock myScopedLock(mut);
    //StartCriticalSection();
    int unusedID = uniqueIDgenerator->generateUniqueID(); //generateUnitID();
    BoxUnit unit(B, unusedID,generateLocalID());
    boxUnits.push_back(unit);
    setSelectedUnitAndBox(unusedID, 0);
    //EndCriticalSection();
    return unusedID;
}

void SpikeSortBoxes::getUnitColor(int UnitID, uint8& R, uint8& G, uint8& B)
{
    for (int k = 0; k < boxUnits.size(); k++)
    {
        if (boxUnits[k].getUnitID() == UnitID)
        {
            R = boxUnits[k].ColorRGB[0];
            G = boxUnits[k].ColorRGB[1];
            B = boxUnits[k].ColorRGB[2];
            break;
        }
    }
    std::vector<PCAUnit> &pcaUnits = pca_sorting_.units;
    for (int k = 0; k < pcaUnits.size(); k++)
    {
        if (pcaUnits[k].getUnitID() == UnitID)
        {
            R = pcaUnits[k].ColorRGB[0];
            G = pcaUnits[k].ColorRGB[1];
            B = pcaUnits[k].ColorRGB[2];
            break;
        }
    }
}

int SpikeSortBoxes::generateLocalID()
{
    // finds the first unused ID and return it

    int ID = 1;

    while (true)
    {
        bool used=false;
        for (int k = 0; k < boxUnits.size(); k++)
        {
            if (boxUnits[k].getLocalID() == ID)
            {
                used = true;
                break;
            }
        }
        std::vector<PCAUnit> &pcaUnits = pca_sorting_.units;
        for (int k = 0; k < pcaUnits.size(); k++)
        {
            if (pcaUnits[k].getLocalID() == ID)
            {
                used = true;
                break;
            }
        }

        if (used)
            ID++;
        else
            break;
    }
    return ID;
}

int SpikeSortBoxes::generateUnitID()
{

    int ID = uniqueIDgenerator->generateUniqueID();
    return ID;
}


void SpikeSortBoxes::generateNewIDs()
{
    const ScopedLock myScopedLock(mut);
    for (int k=0; k<boxUnits.size(); k++)
    {
        boxUnits[k].UnitID = generateUnitID();
    }
    for (int k=0; k<pca_sorting_.units.size(); k++)
    {
        pca_sorting_.units[k].UnitID = generateUnitID();
    }
    synchronizePCAParameter();
}

void SpikeSortBoxes::removeAllUnits()
{
    const ScopedLock myScopedLock(mut);
    boxUnits.clear();
    pca_sorting_.units.clear();
    synchronizePCAParameter();
}

bool SpikeSortBoxes::removeUnit(int unitID)
{
    const ScopedLock myScopedLock(mut);
    for (int k=0; k<boxUnits.size(); k++)
    {
        if (boxUnits[k].getUnitID() == unitID)
        {
            boxUnits.erase(boxUnits.begin()+k);
            return true;
        }
    }

    std::vector<PCAUnit> &pcaUnits = pca_sorting_.units;
    for (int k=0; k < pcaUnits.size(); k++)
    {
        if (pcaUnits[k].getUnitID() == unitID)
        {
            pcaUnits.erase(pcaUnits.begin() + k);
            synchronizePCAParameter();
            return true;
        }
    }

    return false;

}

bool SpikeSortBoxes::addBoxToUnit(int channel, int unitID)
{
    const ScopedLock myScopedLock(mut);

    //StartCriticalSection();

    for (int k = 0; k < boxUnits.size(); k++)
    {
        if (boxUnits[k].getUnitID() == unitID)
        {
            Box B = boxUnits[k].lstBoxes[boxUnits[k].lstBoxes.size() - 1];
            B.x += 100;
            B.y -= 30;
            B.channel = channel;
            boxUnits[k].addBox(B);
            setSelectedUnitAndBox(unitID, (int) boxUnits[k].lstBoxes.size() - 1);
            // EndCriticalSection();
            return true;
        }
    }
    // EndCriticalSection();
    return false;
}


bool SpikeSortBoxes::addBoxToUnit(int channel, int unitID, Box B)
{
    const ScopedLock myScopedLock(mut);
    //StartCriticalSection();
    for (int k=0; k<boxUnits.size(); k++)
    {
        if (boxUnits[k].getUnitID() == unitID)
        {
            boxUnits[k].addBox(B);
            // EndCriticalSection();
            return true;
        }
    }
    //EndCriticalSection();
    return false;
}

std::vector<BoxUnit> SpikeSortBoxes::getBoxUnits()
{
    //StartCriticalSection();
    const ScopedLock myScopedLock(mut);
    std::vector<BoxUnit> unitsCopy = boxUnits;
    //EndCriticalSection();
    return unitsCopy;
}


std::vector<PCAUnit> SpikeSortBoxes::getPCAUnits()
{
    const ScopedLock myScopedLock(mut);
    return pca_sorting_.units;
}

void SpikeSortBoxes::updatePCAUnits(std::vector<PCAUnit> _units)
{
    const ScopedLock myScopedLock(mut);
    pca_sorting_.units = _units;
    synchronizePCAParameter();
}

void SpikeSortBoxes::updateBoxUnits(std::vector<BoxUnit> _units)
{
    const ScopedLock myScopedLock(mut);
    boxUnits = _units;
}




// tests whether a candidate spike belongs to one of the defined units
bool SpikeSortBoxes::sortSpike(SorterSpikePtr so, bool PCAfirst)
{
    const ScopedLock myScopedLock(mut);

    std::vector<PCAUnit> &pcaUnits = pca_sorting_.units;
    if (PCAfirst)
    {
        for (int k=0; k < pcaUnits.size(); k++)
        {
            if (pcaUnits[k].isWaveFormInsidePolygon(so))
            {
                so->sortedId = pcaUnits[k].getUnitID();
                so->color[0] = pcaUnits[k].ColorRGB[0];
                so->color[1] = pcaUnits[k].ColorRGB[1];
                so->color[2] = pcaUnits[k].ColorRGB[2];
                return true;
            }
        }

        for (int k=0; k<boxUnits.size(); k++)
        {
            if (boxUnits[k].isWaveFormInsideAllBoxes(so))
            {
                so->sortedId = boxUnits[k].getUnitID();
                so->color[0] = boxUnits[k].ColorRGB[0];
                so->color[1] = boxUnits[k].ColorRGB[1];
                so->color[2] = boxUnits[k].ColorRGB[2];
                boxUnits[k].updateWaveform(so);
                return true;
            }
        }
    }
    else
    {

        for (int k=0; k<boxUnits.size(); k++)
        {
            if (boxUnits[k].isWaveFormInsideAllBoxes(so))
            {
                so->sortedId = boxUnits[k].getUnitID();
                so->color[0] = boxUnits[k].ColorRGB[0];
                so->color[1] = boxUnits[k].ColorRGB[1];
                so->color[2] = boxUnits[k].ColorRGB[2];
                boxUnits[k].updateWaveform(so);
                return true;
            }
        }
        for (int k=0; k<pcaUnits.size(); k++)
        {
            if (pcaUnits[k].isWaveFormInsidePolygon(so))
            {
                so->sortedId = pcaUnits[k].getUnitID();
                so->color[0] = pcaUnits[k].ColorRGB[0];
                so->color[1] = pcaUnits[k].ColorRGB[1];
                so->color[2] = pcaUnits[k].ColorRGB[2];
                pcaUnits[k].updateWaveform(so);
                return true;
            }
        }

    }

    return false;
}


bool  SpikeSortBoxes::removeBoxFromUnit(int unitID, int boxIndex)
{
    const ScopedLock myScopedLock(mut);
    //StartCriticalSection();
    for (int k=0; k<boxUnits.size(); k++)
    {
        if (boxUnits[k].getUnitID() == unitID)
        {
            bool s= boxUnits[k].deleteBox(boxIndex);
            setSelectedUnitAndBox(-1,-1);
            //EndCriticalSection();
            return s;
        }

    }

    //EndCriticalSection();
    return false;
}

std::vector<Box> SpikeSortBoxes::getUnitBoxes(int unitID)
{
    std::vector<Box> boxes;
    const ScopedLock myScopedLock(mut);
    //StartCriticalSection();
    for (int k=0; k< boxUnits.size(); k++)
    {
        if (boxUnits[k].getUnitID() == unitID)
        {

            boxes = boxUnits[k].getBoxes();
            // EndCriticalSection();
            return boxes;
        }
    }
    //EndCriticalSection();
    return boxes;
}


int SpikeSortBoxes::getNumBoxes(int unitID)
{
    const ScopedLock myScopedLock(mut);
    // StartCriticalSection();
    for (int k=0; k< boxUnits.size(); k++)
    {
        if (boxUnits[k].getUnitID() == unitID)
        {

            int n =boxUnits[k].getNumBoxes();
            // EndCriticalSection();
            return n;
        }
    }
    //EndCriticalSection();
    return -1;
}

/*
        public void DrawBoxes(BufferedGraphics buf, int channel, int SelectedUnit, int SelectedBox,
            double XScale, double XOffset, double YScale, double YOffset)
        {
            mut.WaitOne();
            int unitcounter=0;
            foreach (BoxUnit unit in lstChannelUnits[channel])
            {
                Pen myPen = new Pen(unit.ColorRGB);
                int boxcounter =0;
                foreach (Box b in unit.getBoxes())
                {
                    if (unit.GetUnitID() == SelectedUnit && SelectedBox == boxcounter)
                    {
                        myPen.Width = 3;

                    }
                    else
                    {
                        myPen.Width = 1;
                    }

                   buf.Graphics.DrawRectangle(myPen, (float)(b.x * XScale + XOffset),
                                                      (float)(YOffset-YScale*b.y),
                                                      (float)((b.w) * XScale),
                                                      (float)(YScale*b.h));
                    boxcounter++;
                }
                unitcounter++;
            }
            mut.ReleaseMutex();
        }

        public bool IsUnitActivated(int channel, int unitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
                if (unit.GetUnitID() == unitID)
                    return unit.IsActivated();
            return false;
        }

        public void ActivateUnit(int channel, int unitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
                if (unit.GetUnitID() == unitID)
                     unit.ActivateUnit();
        }
        public void DectivateUnit(int channel, int unitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
                if (unit.GetUnitID() == unitID)
                    unit.DeactivateUnit();
        }
        public double GetNumSecUnitActive(int channel, int unitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
                if (unit.GetUnitID() == unitID)
                    return unit.GetNumSecondsActive();

            return 0;
        }
        public void ToggleActivate(int channel, int unitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
                if (unit.GetUnitID() == unitID)
                    unit.ToggleActive();
        }


        public void UpdateUnitBoxPos(int channel, int unitID, int boxID, PointD P)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
            if (unit.GetUnitID() == unitID)
                unit.SetBoxPos(boxID,P);
        }

        public void UpdateUnitBoxPosAndSize(int channel, int unitID, int boxID, Box B)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
                if (unit.GetUnitID() == unitID)
                    unit.SetBox(boxID, B);
        }


        //YScaleFactors[yScaleIndex]
        public void DrawWavesInBuffer(BufferedGraphics buf, int channel, double XScale, double YScale, double YOffset)
        {
            mut.WaitOne();
            foreach (WaveForm wf in WaveFormBuffer[channel])
            {
                Pen myPen = new Pen(wf.colorRGB);
                for (int j = 0; j < wf.numpts - 1; j++)
                {
                    buf.Graphics.DrawLine(myPen,
                        (float)(XScale * j),
                        (float)(YOffset - YScale * wf.Wave[wf.detected_on_channel_index, j]),
                        (float)(XScale * (j + 1)), (float)(YOffset - YScale * wf.Wave[wf.detected_on_channel_index, j + 1]));
                }
            }
            mut.ReleaseMutex();
        }

         public void ClearAvgWaveForm(int channel, int UnitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
            {
                if (unit.GetUnitID() == UnitID)
                {
                    unit.ResetWaveForm();
                    return;
                }
            }
        }

        public void GetWaveFormsMeanVar(int Channel, out double[,] Mean, out double[,] Std, out int[] UnitIDs, out Boolean[] HasData)
        {
            mut.WaitOne();
            int Counter=0,dim1;
            if (WaveFormBuffer[0].Count == 0)
                dim1 = 62;
            else
                dim1 = WaveFormBuffer[0][0].Wave.GetLength(1);

            HasData = new Boolean[lstChannelUnits[Channel].Count];
            Mean = new double[lstChannelUnits[Channel].Count,dim1];
            Std =  new double[lstChannelUnits[Channel].Count,dim1];
            UnitIDs = new int[lstChannelUnits[Channel].Count];
            foreach (BoxUnit unit in lstChannelUnits[Channel])
            {
                double[] MeanUnit, StdUnit;
                HasData[Counter] = false;
                UnitIDs[Counter] = unit.GetUnitID();
                unit.GetWaveFormMeanStd(out MeanUnit, out StdUnit, 0);
                if (MeanUnit == null)
                {
                    Counter++;
                    continue;
                }
                HasData[Counter] = true;
                for (int j = 0; j < dim1; j++)
                {
                    Mean[Counter,j] = MeanUnit[j];
                    Std[Counter, j] = StdUnit[j];
                }

                Counter++;
            }
            mut.ReleaseMutex();
        }


        private bool GetUnit(int Channel, int UnitID, out BoxUnit box)
        {
            box = new BoxUnit(0);
            foreach (BoxUnit unit in lstChannelUnits[Channel])
            {
                if (unit.GetUnitID() == UnitID)
                {
                    box = unit;
                    return true;
                }
            }
            return false;
        }

        public void DrawInterSpikeHistogram(int channel, int UnitID, BufferedGraphics buf, int Width, int Height)
        {
            mut.WaitOne();
            BoxUnit unit;
            if (!GetUnit(channel, UnitID, out unit))
            {
                mut.ReleaseMutex();
                return;
            }

            Color FadedColor = Color.FromArgb(unit.ColorRGB.A, unit.ColorRGB.R / 2, unit.ColorRGB.G / 2, unit.ColorRGB.B / 2);
            Pen myPen = new Pen(FadedColor);

            int NumPointsToDraw = 30;
            float XScale = (float)Width / NumPointsToDraw;
            Histogram hist = unit.GetInterSpikeHistogram();
            SolidBrush myBrush = new SolidBrush(FadedColor);
            SolidBrush myRedBrush = new SolidBrush(Color.Red);
            if (hist.Max > 0)
            {
                float YScale = (float)Height / hist.Max;

                 for (int j = 0; j < NumPointsToDraw - 1; j++)
                {
                     if (j > 2)
                        buf.Graphics.FillRectangle(myBrush, new Rectangle((int)(XScale * (j + 0.5)), (int)(Height - YScale * hist.Counter[j]), (int)(Width / NumPointsToDraw), (int)(YScale * hist.Counter[j])));
                     else
                         buf.Graphics.FillRectangle(myRedBrush, new Rectangle((int)(XScale * (j + 0.5)), (int)(Height - YScale * hist.Counter[j]), (int)(Width / NumPointsToDraw), (int)(YScale * hist.Counter[j])));

                }
            }
            mut.ReleaseMutex();
        }


        public void DrawAvgWaveForm(int channel, int UnitID, BufferedGraphics buf, int Width, int Height, double YScale )
        {
            BoxUnit unit;
            mut.WaitOne();
            if (!GetUnit(channel, UnitID, out unit))
            {
                mut.ReleaseMutex();
                return;
            }

            double [] Mean;
            double [] Std;
            unit.GetWaveFormMeanStd(out Mean, out Std, 0);
            if (Mean == null)
            {
                mut.ReleaseMutex();
                return;
            }


            Pen myPen = new Pen(unit.ColorRGB);
            Pen myDashPen = new Pen(unit.ColorRGB);
            myDashPen.DashStyle = System.Drawing.Drawing2D.DashStyle.DashDot;

            int WaveFormLength = Mean.Length;
            float XScale = (float)Width / 62.0F;
            float YOffset = (float)Height / 2.0F;

            for (int j = 0; j < WaveFormLength - 1; j++)
            {
                buf.Graphics.DrawLine(myPen,
                    (float)(XScale * j),
                    (float)(YOffset - YScale * Mean[j]),
                    (float)(XScale * (j + 1)), (float)(YOffset - YScale * Mean[j + 1]));

                buf.Graphics.DrawLine(myDashPen,
                                   (float)(XScale * j),
                                   (float)(YOffset - YScale * (Std[ j] + Mean[ j])),
                                   (float)(XScale * (j + 1)), (float)(YOffset - YScale * (Std[ j + 1] +Mean[ j + 1])));

                buf.Graphics.DrawLine(myDashPen,
                                   (float)(XScale * j),
                                   (float)(YOffset - YScale * (-Std[ j] + Mean[ j])),
                                   (float)(XScale * (j + 1)), (float)(YOffset - YScale * (-Std[j + 1] + Mean[ j + 1])));

            }
            mut.ReleaseMutex();
        }


        public Histogram GetInterSpikeHistogram(int channel, int unitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
            {
                if (unit.GetUnitID() == unitID)
                {
                    return unit.GetInterSpikeHistogram();
                }
            }
            return new Histogram(0,0,0);
        }

        public bool QueryNewData(int channel, int unitID)
        {
            foreach (BoxUnit unit in lstChannelUnits[channel])
            {
                if (unit.GetUnitID() == unitID)
                {
                    return unit.QueryNewData();
                }
            }
            return false;
        }

        public void ClearBuffer(int ch)
        {
            WaveFormBuffer[ch].Clear();
        }

        public WaveForm GetLastWaveForm(int channel)
        {
            mut.WaitOne();
            WaveForm wf ;
            if (WaveFormBuffer[channel].Count > 0)
                wf = WaveFormBuffer[channel][WaveFormBuffer[channel].Count - 1];
            else
                wf = new WaveForm();

            mut.ReleaseMutex();
            return wf;
        }

		*/

/**************************/


bool cPolygon::isPointInside(PointD p) const
{
    PointD p1, p2;

    bool inside = false;

    if (pts_.size() < 3)
    {
        return inside;
    }

    float offset_x = offset_[0];
    float offset_y = offset_[1];

    PointD oldPoint(pts_[pts_.size() - 1][0] + offset_x, pts_[pts_.size() - 1][1] + offset_y);

    for (int i = 0; i < pts_.size(); i++)
    {
        PointD newPoint(pts_[i][0] + offset_x, pts_[i][1] + offset_y);

        if (newPoint[0] > oldPoint[0])
        {
            p1 = oldPoint;
            p2 = newPoint;
        }
        else
        {
            p1 = newPoint;
            p2 = oldPoint;
        }

        if ((newPoint[0] < p[0]) == (p[0] <= oldPoint[0])
            && ((p[1] - p1[1]) * (p2[0] - p1[0]) < (p2[1] - p1[1]) * (p[0] - p1[0]))) {
            inside = !inside;
        }

        oldPoint = newPoint;
    }

    return inside;
}


json cPolygon::ToJson() const {
    json ret;
    ret["type"] = "POLYGON";
    for (const auto &pt: pts_) {
        json pt_json = pt.dims();
        ret["points"].push_back(pt_json);
    }
    return ret;
}

bool cPolygon::FromJson(const json &value) {
    if (!value.contains("type") || !value.contains("points") || !value["points"].is_array()) {
        return false;
    }
    if (value["type"].get<std::string>() != "POLYGON") {
        return false;
    }
    pts_.clear();
    try {
        for (const auto& pt_json : value["points"]) {
            pts_.emplace_back(pt_json.get<std::vector<float>>());
        }
    } catch (json::exception &) {
        return false;
    }
    return true;
}







/*****************/


/*************************/
PCAUnit::PCAUnit()
{
}

PCAUnit::PCAUnit(int ID, int localID_): UnitID(ID),localID(localID_)
{
    BoxUnit::setDefaultColors(ColorRGB, localID);
};

PCAUnit::~PCAUnit()
{
}

int PCAUnit::getUnitID()
{
    return UnitID;
}
int PCAUnit::getLocalID()
{
    return localID;
}

bool PCAUnit::isPointInsidePolygon(PointD p)
{
    return poly->isPointInside(p);
}

bool PCAUnit::isWaveFormInsidePolygon(SorterSpikePtr so)
{
    // Only support 2D projections with these "hand-drawn" polygons
    return poly->isPointInside(PointD(so->pcProj[0],so->pcProj[1]));
}

void PCAUnit::resizeWaveform(int newlength)
{
}


void PCAUnit::updateWaveform(SorterSpikePtr so)
{
    WaveformStat.update(so);
}

json PCAUnit::ToJson() const {
    json ret;
    ret["id"] = UnitID;
    ret["local_id"] = localID;
    ret["shape"] = poly->ToJson();
    return ret;
}

bool PCAUnit::FromJson(const json &value) {
    try {
        UnitID = value["id"];
        localID = value["local_id"];
        if (value["shape"]["type"] == "POLYGON") {
            poly = std::make_shared<cPolygon>();
        } else if (value["shape"]["type"] == "ELLIPSE") {
            poly = std::make_shared<cEllipse>();
        } else {
            std::cout << "Unhandled shape type?? " << value["shape"]["type"] << std::endl;
            return false;
        }
        return poly->FromJson(value["shape"]);
    } catch (json::exception &e) {
        std::cout << "Could not parse json into PCA Unit. Error: " << e.what() << std::endl;
        return false;
    }
}

/***************************/


/*
  An implementation of SVD from Numerical Recipes in C and Mike Erhdmann's lectures
*/


#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1 = (a),maxarg2 = (b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1 = (a),iminarg2 = (b),(iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

PCAjob::PCAjob(SorterSpikeArray& _spikes,
               PCAResults *results,
               std::atomic<bool> &_reportDone)
        : spikes(_spikes), results_(results), reportDone(_reportDone) {
	SorterSpikePtr spike = spikes[0];
    cov = nullptr;

    dim = spike->getChannel()->getNumChannels()*spike->getChannel()->getTotalSamples();

};

PCAjob::~PCAjob()
{

}



// calculates sqrt( a^2 + b^2 ) with decent precision
float PCAjob::pythag(float a, float b)
{
    float absa,absb;

    absa = fabs(a);
    absb = fabs(b);

    if (absa > absb)
        return (absa * sqrt(1.0 + SQR(absb/absa)));
    else
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

/*
  Modified from Numerical Recipes in C
  Given a matrix a[nRows][nCols], svdcmp() computes its singular value
  decomposition, A = U * W * Vt.  A is replaced by U when svdcmp
  returns.  The diagonal matrix W is output as a vector w[nCols].
  V (not V transpose) is output as the matrix V[nCols][nCols].
*/
int PCAjob::svdcmp(float** a, int nRows, int nCols, float* w, float** v)
{

    int flag, i, its, j, jj, k, l = 0, nm = 0;
    float anorm, c, f, g, h, s, scale, x, y, z, *rv1;

    rv1 = new float[nCols];
    if (rv1 == NULL)
    {
        printf("svdcmp(): Unable to allocate vector\n");
        return (-1);
    }

    g = scale = anorm = 0.0;
    for (i = 0; i < nCols; i++)
    {
        l = i+1;
        rv1[i] = scale*g;
        g = s = scale = 0.0;
        if (i < nRows)
        {
            for (k = i; k < nRows; k++)
            {
                //std::cout << k << " " << i << std::endl;
                scale += fabs(a[k][i]);
            }

            if (scale)
            {
                for (k = i; k < nRows; k++)
                {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s),f);
                h = f * g - s;
                a[i][i] = f - g;

                for (j = l; j < nCols; j++)
                {
                    for (s = 0.0, k = i; k < nRows; k++) s += a[k][i] * a[k][j];
                    f = s / h;
                    for (k = i; k < nRows; k++) a[k][j] += f * a[k][i];
                }

                for (k = i; k < nRows; k++)
                    a[k][i] *= scale;
            } // end if (scale)
        } // end if (i < nRows)
        w[i] = scale * g;
        g = s = scale = 0.0;
        if (i < nRows && i != nCols-1)
        {
            for (k = l; k < nCols; k++) scale += fabs(a[i][k]);
            if (scale)
            {
                for (k = l; k < nCols; k++)
                {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = - SIGN(sqrt(s),f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k=l; k<nCols; k++) rv1[k] = a[i][k] / h;
                for (j=l; j<nRows; j++)
                {
                    for (s=0.0,k=l; k<nCols; k++) s += a[j][k] * a[i][k];
                    for (k=l; k<nCols; k++) a[j][k] += s * rv1[k];
                }
                for (k=l; k<nCols; k++) a[i][k] *= scale;
            }
        }
        anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));


    }

    for (i=nCols-1; i>=0; i--)
    {
        if (i < nCols-1)
        {
            if (g)
            {
                for (j=l; j<nCols; j++)
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for (j=l; j<nCols; j++)
                {
                    for (s=0.0,k=l; k<nCols; k++) s += a[i][k] * v[k][j];
                    for (k=l; k<nCols; k++) v[k][j] += s * v[k][i];
                }
            }
            for (j=l; j<nCols; j++) v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    for (i=IMIN(nRows,nCols) - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        for (j=l; j<nCols; j++) a[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            for (j=l; j<nCols; j++)
            {
                for (s=0.0,k=l; k<nRows; k++) s += a[k][i] * a[k][j];
                f = (s / a[i][i]) * g;
                for (k=i; k<nRows; k++) a[k][j] += f * a[k][i];
            }
            for (j=i; j<nRows; j++) a[j][i] *= g;
        }
        else
            for (j=i; j<nRows; j++) a[j][i] = 0.0;
        ++a[i][i];
    }

    for (k=nCols-1; k>=0; k--)
    {
        for (its=0; its<30; its++)
        {
            flag = 1;
            for (l=k; l>=0; l--)
            {
                nm = l-1;
                if ((fabs(rv1[l]) + anorm) == anorm)
                {
                    flag =  0;
                    break;
                }
                if ((fabs(w[nm]) + anorm) == anorm) break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i=l; i<=k; i++)
                {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((fabs(f) + anorm) == anorm) break;
                    g = w[i];
                    h = pythag(f,g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j=0; j<nRows; j++)
                    {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            z = w[k];
            if (l == k)
            {
                if (z < 0.0)
                {
                    w[k] = -z;
                    for (j=0; j<nCols; j++) v[j][k] = -v[j][k];
                }
                break;
            }
            //if(its == 29) printf("no convergence in 30 svdcmp iterations\n");
            x = w[l];
            nm = k-1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f,1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f))) - h)) / x;
            c = s = 1.0;
            for (j=l; j<=nm; j++)
            {
                i = j+1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f,h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj=0; jj<nCols; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = pythag(f,h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj=0; jj < nRows; jj++)
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    delete[] rv1;

    return (0);
}


void PCAjob::computeCov()
{
    // allocate and zero
    cov = new float*[dim];
    float* mean  = new float[dim];
    for (int k = 0; k < dim; k++)
    {
        cov[k] = new float[dim];
        for (int j=0; j<dim; j++)
        {
            cov[k][j] = 0;
        }
    }
    // compute mean

    for (int j=0; j<dim; j++)
    {
        mean[j] = 0;
        for (int i=0; i<spikes.size(); i++)
        {
            SorterSpikePtr spike = spikes[i];
            float v = spikeDataIndexToMicrovolts(spike, j) ;
            mean[j] += v / dim;
        }
    }
    // aggregate


    for (int i=0; i<dim; i++)
    {
        for (int j=i; j<dim; j++)
        {
            // cov[i][j] = sum_k[ (X(i,:)) * (Xj-mue(j) ]
            float sum = 0 ;
            for (int k=0; k<spikes.size(); k++)
            {

                SorterSpikePtr spike = spikes[k];
                float vi = spikeDataIndexToMicrovolts(spike, i);
                float vj = spikeDataIndexToMicrovolts(spike, j);
                sum += (vi-mean[i]) * (vj-mean[j]);
            }
            cov[i][j] = sum / (dim-1);
            cov[j][i] = sum / (dim-1);
        }
    }
    delete[] mean;

    // delete covariances
    //for (int k = 0; k < dim; k++)
    //delete cov[k];

    //delete(cov);
    //cov = nullptr;

}

std::vector<int> sort_indexes(std::vector<float> v)
{
    // initialize original index locations
    std::vector<int> idx(v.size());

    for (int i = 0; i != idx.size(); ++i)
    {
        idx[i] = i;
    }

    //sort indexes based on comparing values in v
    sort(
        idx.begin(),
        idx.end()//,
        //[&v](size_t i1, size_t i2)
        //{
        //	return v[i1] > v[i2];
        //}
    );

    return idx;
}

void PCAjob::computeSVD()
{


    float** eigvec, *sigvalues;
    sigvalues = new float[dim];

    eigvec = new float*[dim];
    for (int k = 0; k < dim; k++)
    {
        eigvec[k] = new float[dim];
        for (int j=0; j<dim; j++)
        {
            eigvec[k][j] = 0;
        }
    }

    svdcmp(cov, dim, dim, sigvalues, eigvec);

    std::vector<float> sig;
    sig.resize(dim);
    for (int k = 0; k < dim; k++)
        sig[k] = sigvalues[k];

    std::vector<int> sortind = sort_indexes(sig);

    results_->clear();
    int num_pcs = results_->num_pcs;

    for (int pc_idx = 0; pc_idx < num_pcs; pc_idx++) {
        std::vector<float> pc;
        for (int k = 0; k < dim; k++)
        {
            pc.push_back(eigvec[k][sortind[pc_idx]]);
        }
        results_->pcs.push_back(pc);
    }

    // project samples to find the display range
    for (int pc_idx = 0; pc_idx < num_pcs; pc_idx++) {
        float min = 1e10;
        float max = -1e10;
        for (auto spike : spikes) {
            float sum = 0;
            for (int k = 0; k < dim; k++) {
                sum += spikeDataIndexToMicrovolts(spike, k) * results_->pcs[pc_idx][k];
            }
            min = std::min(min, sum);
            max = std::max(max, sum);
        }

        results_->pc_mins.push_back(min - 1.5 * (max - min));
        results_->pc_maxes.push_back(max + 1.5 * (max - min));
    }

    // clear memory
    for (int k = 0; k < dim; k++)
    {
        delete[] eigvec[k];
    }
    delete[] eigvec;
    delete[] sigvalues;

    // delete covariances
    for (int k = 0; k < dim; k++)
        delete[] cov[k];

    delete[] cov;
    cov = nullptr;

}


/**********************/


void PCAcomputingThread::addPCAjob(PCAJobPtr job)
{
	{
		ScopedLock critical(lock);
		jobs.add(job);
	}
	
    if (!isThreadRunning())
    {
        startThread();
    }
}

void PCAcomputingThread::run()
{
    while (jobs.size() > 0)
    {
		lock.enter();
        PCAJobPtr J = jobs.removeAndReturn(0);
	if (J == nullptr) continue;
		lock.exit();
        // compute PCA
        // 1. Compute Covariance matrix
        // 2. Apply SVD on covariance matrix
        // 3. Extract the two principal components corresponding to the largest singular values

        J->computeCov();
        J->computeSVD();

        // 4. Report to the spike sorting electrode that PCA is finished
        J->reportDone = true;
    }
}


PCAcomputingThread::PCAcomputingThread() : Thread("PCA")
{

}


/**************************/

float spikeDataBinToMicrovolts(SorterSpikePtr s, int bin, int ch)
{
	jassert(ch >= 0 && ch < s->getChannel()->getNumChannels());
	jassert(bin >= 0 && bin <= s->getChannel()->getTotalSamples());
	float v = s->getData()[bin + ch*s->getChannel()->getTotalSamples()];
	return v;
}


float spikeDataIndexToMicrovolts(SorterSpikePtr s, int index)
{
	float v = s->getData()[index];
	return v;
}

float spikeTimeBinToMicrosecond(SorterSpikePtr s, int bin, int ch)
{
	float spikeTimeSpan = 1.0f / s->getChannel()->getSampleRate() * s->getChannel()->getTotalSamples() * 1e6;
	return float(bin) / (s->getChannel()->getTotalSamples() - 1) * spikeTimeSpan;
}

int microSecondsToSpikeTimeBin(SorterSpikePtr s, float t, int ch)
{
	// Lets say we have 32 samples per wave form

	// t = 0 corresponds to the left most index.
	float spikeTimeSpan = (1.0f / s->getChannel()->getSampleRate() * s->getChannel()->getTotalSamples())*1e6;
	return MIN(s->getChannel()->getTotalSamples() - 1, MAX(0, t / spikeTimeSpan * (s->getChannel()->getTotalSamples() - 1)));
}


SorterSpikeContainer::SorterSpikeContainer(const SpikeChannel* channel, SpikeEvent::SpikeBuffer& spikedata, int64 timestamp)
{
    color[0] = color[1] = color[2] = 127;

    // Hard-code at least 2 projections to allow drawing of bare axes
    pcProj.clear();
	pcProj.resize(2, 0);
	sortedId = 0;
	this->timestamp = timestamp;
	chan = channel;
	int nSamples = chan->getNumChannels() * chan->getTotalSamples();
	data.malloc(nSamples);
	memcpy(data.getData(), spikedata.getRawPointer(), nSamples*sizeof(float));
}

const float* SorterSpikeContainer::getData() const
{
	return data.getData();
}

const SpikeChannel* SorterSpikeContainer::getChannel() const
{
	return chan;
}

int64 SorterSpikeContainer::getTimestamp() const
{
	return timestamp;
}

json PCAResults::ToJson() const {
    json ret;
    ret["num_pcs"] = num_pcs;
    ret["pc_mins"] = pc_mins;
    ret["pc_maxes"] = pc_maxes;
    for (const auto& pc_vec : pcs) {
        ret["pcs"].push_back(pc_vec);
    }
    return ret;
}

bool PCAResults::FromJson(const json &value, int expected_pc_vector_samples) {
    try {
        num_pcs = value["num_pcs"];
    } catch (json::exception &e) {
        std::cout << "Could not parse num_pcs as integer: " << e.what() << std::endl;
        return false;
    }

    try {
        pc_mins = std::vector<float>(value["pc_mins"].begin(), value["pc_mins"].end());
    } catch (json::exception &e) {
        std::cout << "Could not parse pc_mins as array of floats: " << e.what() << std::endl;
        return false;
    }

    try {
        pc_maxes = std::vector<float>(value["pc_maxes"].begin(), value["pc_maxes"].end());
    } catch (json::exception &e) {
        std::cout << "Could not parse pc_maxes as array of floats: " << e.what() << std::endl;
        return false;
    }

    pcs.clear();
    try {
        for (std::vector<float> pc : value["pcs"]) {
            if (pc.size() != expected_pc_vector_samples) {
                std::cout << "Invalid PC vector length given. Need length " << expected_pc_vector_samples << std::endl;
                return false;
            }
            pcs.push_back(pc);
        }
    } catch (json::exception &e) {
        std::cout << "Could not parse pcs as array of array of floats: " << e.what() << std::endl;
        return false;
    }

    if (num_pcs < 2) {
        std::cout << "Invalid number of PCs given. " << num_pcs << std::endl;
        return false;
    }

    // Only accept the results if it populated these fully!
    if (is_populated()) {
        std::cout << "Parsed new PCA results: " << (std::string) *this << std::endl;
        return true;
    }
    return false;
}

cEllipse::cEllipse(const PointD &center,
                   const std::vector<std::vector<float>> &rotation,
                   const std::vector<float> &radii) :
        center_(center), rotation_(rotation), radii_(radii) {
    jassert(!rotation.empty());
    // Rotation is NxN, radii is an Nx1 vector, and Center is an N-dim point.
    jassert(rotation.size() == radii.size());
    jassert(rotation.size() == center.dims().size());

    n_dims = radii.size();

    populate_sample_pts();
}

void cEllipse::populate_sample_pts() {
    int n_sample_pts = 20;
    // Create N points around a unit circle in 2D, and then transform to make them span an ellipse.
    for (int sample_pt_idx = 0; sample_pt_idx < n_sample_pts; sample_pt_idx++) {
        double rads = 2 * M_PI * sample_pt_idx / (double) n_sample_pts;
        std::vector<float> dims;
        dims.push_back(cos(rads));
        dims.push_back(sin(rads));
        for (int j = 2; j < n_dims; j++) {
            dims.push_back(0);
        }

        PointD point = PointD(dims);

        // Stretch according to the radii
        point = point * radii_;

        std::vector<float> rotated;
        rotated.resize(n_dims, 0.0);
        for (int i = 0; i < n_dims; i++) {
            float res = 0;
            for (int j = 0; j < n_dims; j++) {
                res += rotation_[j][i] * point[j];
            }
            rotated[i] = res;
        }
        // Rotate...
        point.set(rotated);

        // Then translate to center of the ellipse:
        point += center_;

        sample_pts_.push_back(point);
        for (const float dim : point.dims()) {
            std::cout << dim << " ";
        }
        std::cout << std::endl;
    }
}

bool cEllipse::isPointInside(PointD p) const {
    PointD centered = p - center_;
    // || rotation matrix * (pt - center) / radii ||_2 < 1?
    double norm2 = 0;
    for (int i = 0; i < n_dims; i++) {
        float res = 0;
        for (int j = 0; j < n_dims; j++) {
            res += rotation_[i][j] * centered[j];
        }
        norm2 += std::pow(res / radii_[i], 2);
    }

    return norm2 < 1;
}

json cEllipse::ToJson() const {
    json ret;
    ret["type"] = "ELLIPSE";
    ret["center"] = center_.dims();
    ret["radii"] = radii_;

    json rotation_json;
    for (const auto& rotation_row : rotation_) {
        rotation_json.push_back(rotation_row);
    }
    ret["rotation"] = rotation_json;
    return ret;
}

bool cEllipse::FromJson(const json &value) {
    if (!value.contains("type") ||
        value["type"].get<std::string>() != "ELLIPSE" ||
        !value.contains("center") ||
        !value.contains("radii") ||
        !value.contains("rotation")) {
        return false;
    }

    rotation_.clear();
    try {
        for (const auto& rotation_row_json : value["rotation"]) {
            rotation_.push_back(rotation_row_json.get<std::vector<float>>());
        }
    } catch (json::exception &) {
        return false;
    }

    try {
        center_ = PointD(value["center"].get<std::vector<float>>());
        radii_ = value["radii"].get<std::vector<float>>();
    } catch (json::exception &) {
        return false;
    }

    if (center_.dims().size() != rotation_.size() || radii_.size() != rotation_.size()) {
        return false;
    }
    n_dims = rotation_.size();
    populate_sample_pts();
    return true;
}

Value PCASorting::ToValue() const {
    json sorting;
    sorting["results"] = results.ToJson();
    sorting["pca_units"] = json::array();
    for (const auto &unit : units) {
        sorting["pca_units"].push_back(unit.ToJson());
    }
    return Value(var(sorting.dump()));
}

bool PCASorting::FromValue(const Value &value, int expected_pc_vector_samples) {
    juce::String parameter_value_str = value.toString();

    json sorting;
    try {
        sorting = json::parse(parameter_value_str.toStdString());
    } catch (json::exception &e) {
        std::cout << "Expected json in electrode parameter, but failed to parse: " << e.what() << std::endl;
        return false;
    }

    results = PCAResults();
    bool success = results.FromJson(sorting["results"], expected_pc_vector_samples);
    if (!success) {
        return false;
    }

    units.clear();
    for (const auto &unit_json : sorting["pca_units"]) {
        PCAUnit unit;
        success = unit.FromJson(unit_json);
        if (!success) {
            return false;
        }
        units.push_back(unit);
    }
    return true;
}
