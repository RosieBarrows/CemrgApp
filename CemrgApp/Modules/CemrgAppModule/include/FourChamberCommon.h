#ifndef FOURCHAMBERCOMMON_H
#define FOURCHAMBERCOMMON_H

#include <QString>
#include <QStringList>
#include <QJsonObject>
#include <QJsonArray>

enum ManualPoints { CYLINDERS, SLICERS, VALVE_PLAINS };
enum LabelsType {  
    BACKGROUND = 0,
    BLOODPOOL = 1, 
    LEFT_VENTRICLE = 2, 
    RIGHT_VENTRICLE = 3, 
    LEFT_ATRIUM = 4, 
    RIGHT_ATRIUM = 5, 
    AORTA = 6,
    PULMONARY_ARTERY = 7, 
    LSPV = 8, 
    LIPV = 9, 
    RSPV = 10, 
    RIPV = 11, 
    LAA = 12,
    SVC = 13,
    IVC = 14
};
class SegmentationLabels {
    private:
        std::unordered_map<LabelsType, unsigned int> labelMap;

        void UpdateLabels(const SegmentationLabels& other) {
            for (const auto& pair : other.labelMap) {
                labelMap[pair.first] = pair.second; // Update the labelMap
            }
        }


    public:
        SegmentationLabels() {
            labelMap = {
                {BACKGROUND, 0},
                {BLOODPOOL, 1},
                {LEFT_VENTRICLE, 2},
                {RIGHT_VENTRICLE, 3},
                {LEFT_ATRIUM, 4},
                {RIGHT_ATRIUM, 5},
                {AORTA, 6},
                {PULMONARY_ARTERY, 7},
                {LSPV, 8},
                {LIPV, 9},
                {RSPV, 10},
                {RIPV, 11},
                {LAA, 12},
                {SVC, 13},
                {IVC, 14}};
        }

        SegmentationLabels(const SegmentationLabels& other) {
            Update(other);
        }

        unsigned int Get(LabelsType labelType) const {
            auto it = labelMap.find(labelType);
            if (it != labelMap.end()) {
                return it->second;
            }
            return 0; // Default to BACKGROUND
        }

        unsigned int GetLabelFromString(const std::string &labelNameString) const {
            for (const auto &pair : labelMap) {
                if (LabelName(pair.first) == labelNameString) {
                    return pair.second;
                }
            }
            return 0; // Default to BACKGROUND
        }

        void Set(LabelsType labelType, unsigned int value) {
            labelMap[labelType] = value;
        }

        void SetLabelFromString(const std::string &labelNameString, unsigned int newTag) {
            for (auto &pair : labelMap) {
                if (LabelName(pair.first) == labelNameString) {
                    pair.second = newTag;
                    return; // Exit the loop once the label is found and updated
                }
            }
        }

        void Update(const SegmentationLabels& other) {
            std::unordered_map<LabelsType, unsigned int> otherMap;
            other.GetMap(otherMap); // Retrieve the labelMap from the other instance
            labelMap = otherMap;    // Synchronize the maps
        }

        void GetMap(std::unordered_map<LabelsType, unsigned int>& map) const {
            map = labelMap; // Fill the provided map with the unsigned internal labelMap
        }

        std::string LabelName(LabelsType labelType) const {
            switch (labelType) {
            case BACKGROUND: return "BACKGROUND";
            case BLOODPOOL: return "BLOODPOOL";
            case LEFT_VENTRICLE: return "LEFT_VENTRICLE";
            case RIGHT_VENTRICLE: return "RIGHT_VENTRICLE";
            case LEFT_ATRIUM: return "LEFT_ATRIUM";
            case RIGHT_ATRIUM: return "RIGHT_ATRIUM";
            case AORTA: return "AORTA";
            case PULMONARY_ARTERY: return "PULMONARY_ARTERY";
            case LSPV: return "LSPV";
            case LIPV: return "LIPV";
            case RSPV: return "RSPV";
            case RIPV: return "RIPV";
            case LAA: return "LAA";
            case SVC: return "SVC";
            case IVC: return "IVC";
            default: return "UNKNOWN";
            }
        }


        bool LabelExists(unsigned int label) const {
            for (const auto &pair : labelMap) {
                if (pair.second == label) {
                    return true;
                }
            }
            return false;
        }

        unsigned int GetLargestLabel() const {
            unsigned int largestLabel = 0;
            for (const auto &pair : labelMap) {
                if (pair.second > largestLabel) {
                    largestLabel = pair.second;
                }
            }
            return largestLabel;
        }

        unsigned int GenerateNewLabel() const {
            return GetLargestLabel() + 1;
        }

        void LoadJsonObject(QJsonObject json) {
            QStringList keys = json.keys();
            MITK_INFO << "Loading json object.";
            for (const auto &key : keys) {
                unsigned int value = json[key].toInt();
                std::cout << "Key: " << key.toStdString() << " Value: " << value << '\n';
                SetLabelFromString(key.toStdString(), json[key].toInt());
            }
        }

        // functions to convert to json
        void LabelInfoLists(QStringList& labelsKeys, QStringList& labelsValues, QStringList& labelsTypes) {
            for (const auto &pair : labelMap) {
                labelsKeys << QString::fromStdString(LabelName(pair.first));
                labelsValues << QString::number(pair.second);
                labelsTypes << "int";
            }
        }

        QStringList LabelNames() const {
            QStringList labelNames;
            for (const auto &pair : labelMap) {
                labelNames << QString::fromStdString(LabelName(pair.first));
            }
            return labelNames;
        }

        QStringList LabelTags() const {
            QStringList labelTags;
            for (const auto &pair : labelMap) {
                labelTags << QString::number(pair.second);
            }
            return labelTags;
        }

        SegmentationLabels& operator=(const SegmentationLabels& other) {
            if (this != &other) {
                Update(other);
            }
            return *this;
        }
};

struct FourChamberSubfolders {
    QString SEG, MESH, UVC, UVC_LA, UVC_RA, AFIB, PRESIM, SIMS, PAR;
    FourChamberSubfolders()
    :SEG("segmentations"),
     MESH("meshing"),
     UVC("surfaces_uvc"),
     UVC_LA("surfaces_uvc_LA"),
     UVC_RA("surfaces_uvc_RA"),
     AFIB("atrial_fibres"),
     PRESIM("pre_simulation"),
     SIMS("sims_folder"), 
     PAR("parfiles")
    {}

    QStringList Subdirectories(){
         return (QStringList() << SEG << MESH << UVC << UVC_LA << UVC_RA << AFIB << PRESIM << SIMS << PAR);
    }
};
struct FourChamberSegmentationNames {
    QString _base, _s1, _s2a, _s2b, _s2c, _s2d, _s2e, _s2f;
    FourChamberSegmentationNames()
        :_base(""), 
         _s1("seg_corrected"), 
         _s2a("seg_s2a"),
         _s2b("seg_s2b"),
         _s2c("seg_s2c"),
         _s2d("seg_s2d"),
         _s2e("seg_s2e"),
         _s2f("seg_s2f") {}
    
    void SetBase(QString name) { _base = name; }
    
    QString QGetBase(QString ext="") { return _base + ext; }
    QString QGetS1(QString ext="") { return _s1 + ext; }
    QString QGetS2A(QString ext="") { return _s2a + ext; }
    QString QGetS2B(QString ext="") { return _s2b + ext; }
    QString QGetS2C(QString ext="") { return _s2c + ext; }
    QString QGetS2D(QString ext="") { return _s2d + ext; }
    QString QGetS2E(QString ext="") { return _s2e + ext; }
    QString QGetS2F(QString ext="") { return _s2f + ext; }

    QString QbaseNii() { return QGetBase(".nii"); }
    QString Qs1Nii() { return QGetS1(".nii"); }
    QString Qs2aNii() { return QGetS2A(".nii"); } 
    QString Qs2bNii() { return QGetS2B(".nii"); } 
    QString Qs2cNii() { return QGetS2C(".nii"); } 
    QString Qs2dNii() { return QGetS2D(".nii"); } 
    QString Qs2eNii() { return QGetS2E(".nii"); } 
    QString Qs2fNii() { return QGetS2F(".nii"); } 

    std::string base() { return QGetBase().toStdString(); }
    std::string s1() { return QGetS1().toStdString(); }
    std::string s2a() { return QGetS2A().toStdString(); } 
    std::string s2b() { return QGetS2B().toStdString(); } 
    std::string s2c() { return QGetS2C().toStdString(); } 
    std::string s2d() { return QGetS2D().toStdString(); } 
    std::string s2e() { return QGetS2E().toStdString(); } 
    std::string s2f() { return QGetS2F().toStdString(); } 

    std::string base_nii() { return QbaseNii().toStdString(); }
    std::string s1_nii() { return Qs1Nii().toStdString(); } 
    std::string s2a_nii() { return Qs2aNii().toStdString(); }
    std::string s2b_nii() { return Qs2bNii().toStdString(); }
    std::string s2c_nii() { return Qs2cNii().toStdString(); }
    std::string s2d_nii() { return Qs2dNii().toStdString(); }
    std::string s2e_nii() { return Qs2eNii().toStdString(); }
    std::string s2f_nii() { return Qs2fNii().toStdString(); }

};

struct ManualPointsStruct {
    QStringList cylinders, slicers, valve_plains;
    ManualPointsStruct(){
        cylinders << "SVC_1" << "SVC_2" << "SVC_3" << "IVC_1" << "IVC_2" << "IVC_3"  
                 << "Ao_1"  << "Ao_2"  << "Ao_3" << "PArt_1" << "PArt_2" << "PArt_3";
        slicers << "SVC_slicer_1" << "SVC_slicer_2" << "SVC_slicer_3" << "SVC_tip"
                << "IVC_slicer_1" << "IVC_slicer_2" << "IVC_slicer_3" << "IVC_tip"
                << "Ao_tip" << "PArt_tip";
        valve_plains << "Ao_WT_tip" << "PArt_WT_tip";
    }

    QStringList CYLINDERS(){ return cylinders; };
    QStringList SLICERS(){ return slicers; };
    QStringList VALVE_PLAINS() { return valve_plains; };

    QStringList GetPointLabelOptions(ManualPoints mpt){
        QStringList res = QStringList();

        switch (mpt) {
            case ManualPoints::CYLINDERS :
                res = CYLINDERS();
                break;

            case ManualPoints::SLICERS : 
                res = SLICERS();
                break;

            case ManualPoints::VALVE_PLAINS :
                res = VALVE_PLAINS();
                break;
            default:
                break;
            return res;
        }
        return res;
    }

    QString title(ManualPoints mpt) {
        QString res = QString();
        switch (mpt) {
            case ManualPoints::CYLINDERS :
                res = "Cylinders";
                break;

            case ManualPoints::SLICERS : 
                res = "Slicers";
                break;

            case ManualPoints::VALVE_PLAINS :
                res = "Valve Plains";
                break;
            default:
                break;
            return res;
        }
        return res;
    }
};

/// @brief Meshtool3D static libraries parameters structures
struct M3DParameters {
    QString seg_dir, seg_name, out_dir, out_name;

    bool mesh_from_segmentation, boundary_relabelling, verbose, eval_thickness;

    double facet_angle, facet_size, facet_distance;
    double cell_rad_edge_ratio, cell_size;
    double rescaleFactor, abs_tol, rel_tol;
    int itr_max, dimKrilovSp;
    
    bool  out_medit, out_carp, out_carp_binary, out_vtk, out_vtk_binary, out_potential;

    M3DParameters()
    :seg_dir(""),
     seg_name(""),
     out_dir(""),
     out_name(""),
     mesh_from_segmentation(true),
     boundary_relabelling(false),
     verbose(true),
     eval_thickness(false),
     facet_angle(30.0),
     facet_size(0.8),
     facet_distance(4),
     cell_rad_edge_ratio(2.0),
     cell_size(0.8),
     rescaleFactor(1000),
     abs_tol(1e-6),
     rel_tol(1e-6),
     itr_max(700),
     dimKrilovSp(500),
     out_medit(false),
     out_carp(true),
     out_carp_binary(false),
     out_vtk(false),
     out_vtk_binary(false),
     out_potential(false)
    {}

    bool CheckFormats(){
        return !(out_medit || out_carp || out_carp_binary ||out_vtk || out_vtk_binary);
    }

    void SetBinaryFormats(){
        out_carp_binary = out_carp;
        out_vtk_binary = out_vtk;
    }

    void KeysAndValues(QStringList& keys, QStringList& values, QStringList& types){
        keys << "seg_dir" << "seg_name" << "mesh_from_segmentation" << "boundary_relabelling" << 
                "facet_angle" << "facet_size" << "facet_distance" << "cell_rad_edge_ratio" << "cell_size" << "rescaleFactor" << 
                "abs_tol" << "rel_tol" << "itr_max" << "dimKrilovSp" << "verbose" << "eval_thickness" << 
                "out_dir" << "out_name" << "out_medit" << "out_carp" << "out_carp_binary" << "out_vtk" << "out_vtk_binary" << "out_potential";
        values << seg_dir << seg_name << QString::number(mesh_from_segmentation) << QString::number(boundary_relabelling)
               << QString::number(facet_angle) << QString::number(facet_size) << QString::number(facet_distance) << QString::number(cell_rad_edge_ratio) 
               << QString::number(cell_size) << QString::number(rescaleFactor) << QString::number(abs_tol)
               << QString::number(rel_tol) << QString::number(itr_max) << QString::number(dimKrilovSp) << QString::number(verbose) << QString::number(eval_thickness) 
               << out_dir << out_name
               << QString::number(out_medit) << QString::number(out_carp) << QString::number(out_carp_binary) << QString::number(out_vtk) << QString::number(out_vtk_binary) << QString::number(out_potential);
        types << "string" << "string" << "int" << "int" << 
                "float" << "float" << "float" << "float" << "float" << "int" << 
                "float" << "float" << "int" << "int" << "int" << "int" << 
                "string" << "string" << "int" << "int" << "int" << "int" << "int" << "int";
    }

};

struct VFibresParams {
    double _alpha_endo, _alpha_epi, _beta_endo, _beta_epi;
    QString _bad_elem, _type, _apex_to_base, _epi, _lv, _rv;
    QString _directory, _meshname, _output;
    int _nonmyo;
     _nomyo_set;

    VFibresParams()
    :_alpha_endo(60),
     _alpha_epi(-60),
     _beta_endo(-65), 
     _beta_epi(25), 
     _bad_elem(""), 
     _type("biv"), 
     _apex_to_base(""), 
     _epi(""),
     _lv(""),
     _rv(""), 
     _nonmyo(-1), 
     _nomyo_set(false)
     {}

    inline void SetDirectory(QString dir) { _directory = dir; };
    inline void SetMeshname(QString name) { _meshname = name; };
    inline void SetOutput(QString out) { _output = out; };

    inline void SetBadElem(QString bad) { _bad_elem = bad; };
    inline void SetApexToBase(QString apex) { _apex_to_base = apex; };

    inline void SetEpi(QString epi) { _epi = epi; };
    inline void SetLV(QString lv) { _lv = lv; };
    inline void SetRV(QString rv) { _rv = rv; };

    inline void SetNonmyo(int nonmyo) { _nonmyo = nonmyo; _nomyo_set = true; };
    inline void SetType(QString type) { _type = type; };

    inline void SetAlphaEndo(double alpha) { _alpha_endo = alpha; };
    inline void SetAlphaEpi(double alpha) { _alpha_epi = alpha; };
    inline void SetBetaEndo(double beta) { _beta_endo = beta; };
    inline void SetBetaEpi(double beta) { _beta_epi = beta; };

    inline QString a_endo() { return QString::number(_alpha_endo); }; 
    inline QString a_epi() { return QString::number(_alpha_epi); }; 
    inline QString b_endo() { return QString::number(_beta_endo); }; 
    inline QString b_epi() { return QString::number(_beta_epi); };
    inline QString nonmyo() { return QString::number(_nonmyo); };
    inline QString type() { return _type; };
    inline QString output() { return directory(_output); };
    inline QString bad_elem(QString dir="") { return directory(_bad_elem); };
    inline QString apex_to_base(QString dir="") { return directory(_apex_to_base); };
    inline QString meshname(QString dir="") { return directory(_meshname); };
    inline QString epi(QString dir="") { return directory(_epi); };
    inline QString lv(QString dir="") { return directory(_lv); };
    inline QString rv(QString dir="") { return directory(_rv); };

    inline QString directory(QString fname = "") { return _directory + "/" + fname; };  

    void KeysAndValues(QStringList& keys, QStringList& values, QStringList& types) {
        keys << "a_endo" << "a_epi" << "b_endo" << "b_epi" << "bad_elem" << "type" << "apex_to_base" 
            << "epi" << "lv" << "rv" << "nonmyo" << "directory" << "meshname" << "output";
        values << a_endo() << a_epi() << b_endo() << b_epi() << bad_elem() << type() << apex_to_base() 
            << epi() << lv() << rv() << nonmyo() << directory() << meshname() << output();
        types << "float" << "float" << "float" << "float" << "string" << "string" << "string" 
            << "string" << "string" << "string" << "int" << "string" << "string" << "string";
    }
};

enum PointsNamesType {
    SVC1, SVC2, SVC3, IVC1, IVC2, IVC3, Ao1, Ao2, Ao3, PArt1, PArt2, PArt3,
    SVC_SLICER1, SVC_SLICER2, SVC_SLICER3, SVC_TIP, IVC_SLICER1, IVC_SLICER2, IVC_SLICER3, IVC_TIP, Ao_TIP, PArt_TIP,
    Ao_WT_TIP, PArt_WT_TIP
};

class BasePointsType {

public:
    BasePointsType(int numEnumValues)
        : points(3 * numEnumValues), pointsSet(false) {}
    
    BasePointsType(ManualPoints mpt)
        : points(3 * GetNumEnumValues(mpt)), pointsSet(false) {}

    double GetPointAt(PointsNamesType ptType, unsigned int index) {
        return points.at(static_cast<int>(ptType) * 3 + index);
    }

    void SetPointAt(PointsNamesType ptType, unsigned int index, double value) {
        points.at(static_cast<int>(ptType) * 3 + index) = value;
    }

    void SetPoint(PointsNamesType ptType, std::vector<double> value) {
        for (int ix = 0; ix < 3; ix++) {
            SetPointAt(ptType, ix, value.at(ix));
        }
    }

    virtual PointsNamesType FromKey(QString key) = 0;

    void SetPoint(QJsonObject json, QString key) {
        if (json[key].isUndefined()) {
            MITK_WARN << ("Undefined key: " + key + '\n').toStdString();
            return;
        }

        for (int ix = 0; ix < 3; ix++) {
            SetPointAt(FromKey(key), ix, json[key].toArray().at(ix).toDouble());
        }
    }

    void SetPointsFromJson(QJsonObject json) {
        QStringList keys = json.keys();
        foreach (QString key, keys) {
            SetPoint(json, key);

        }
    }

    inline void PointSet(bool value) { pointsSet = value; };
    inline void PointSetOn() { PointSet(true); };
    inline void PointSetOff() { PointSet(false); };
    inline bool IsPointSet() { return pointsSet; };

protected:
    std::vector<double> points;
    bool pointsSet;

private: 
    int GetNumEnumValues(ManualPoints mpt) {
        switch (mpt) {
            case ManualPoints::CYLINDERS:
                return static_cast<int>(PArt3) - static_cast<int>(SVC1) + 1; // Or any calculation that fits your needs
            case ManualPoints::SLICERS:
                return static_cast<int>(PArt_TIP) - static_cast<int>(SVC_SLICER1) + 1; // Calculate for SlicersPointsType
            case ManualPoints::VALVE_PLAINS:
                return static_cast<int>(PArt_WT_TIP) - static_cast<int>(Ao_WT_TIP) + 1; // Calculate for ValvePlainsPointsType
        }
    }
};

class CylinderPointsType : public BasePointsType {
public:
    CylinderPointsType()
        : BasePointsType(static_cast<int>(PArt3) - static_cast<int>(SVC1) + 1) {} 

    PointsNamesType FromKey(QString key) override {
        if (key == "SVC_1") return SVC1;
        if (key == "SVC_2") return SVC2;
        if (key == "SVC_3") return SVC3;
        if (key == "IVC_1") return IVC1;
        if (key == "IVC_2") return IVC2;
        if (key == "IVC_3") return IVC3;
        if (key == "Ao_1") return Ao1;
        if (key == "Ao_2") return Ao2;
        if (key == "Ao_3") return Ao3;
        if (key == "PArt_1") return PArt1;
        if (key == "PArt_2") return PArt2;
        if (key == "PArt_3") return PArt3;
        return SVC1; // Default value, you can adjust this as needed
    }
}; 

class SlicersPointsType : public BasePointsType {
public:
    SlicersPointsType()
        : BasePointsType(static_cast<int>(PArt_TIP) - static_cast<int>(SVC_SLICER1) + 1) {} 

    PointsNamesType FromKey(QString key) override {
        if (key == "SVC_slicer_1") return SVC_SLICER1;
        if (key == "SVC_slicer_2") return SVC_SLICER2;
        if (key == "SVC_slicer_3") return SVC_SLICER3;
        if (key == "SVC_tip") return SVC_TIP;
        if (key == "IVC_slicer_1") return IVC_SLICER1;
        if (key == "IVC_slicer_2") return IVC_SLICER2;
        if (key == "IVC_slicer_3") return IVC_SLICER3;
        if (key == "IVC_tip") return IVC_TIP;
        if (key == "Ao_tip") return Ao_TIP;
        if (key == "PArt_tip") return PArt_TIP;
        return SVC_SLICER1; // Default value, you can adjust this as needed
    }
};

class ValvePlainsPointsType : public BasePointsType {
public:
    ValvePlainsPointsType() 
        : BasePointsType(static_cast<int>(PArt_WT_TIP) - static_cast<int>(Ao_WT_TIP) + 1) {} // Total number of enum values for ValvePlains

    PointsNamesType FromKey(QString key) override {
        if (key == "Ao_WT_TIP") return Ao_WT_TIP;
        if (key == "PArt_WT_TIP") return PArt_WT_TIP;
        return Ao_WT_TIP; // Default value, you can adjust this as needed
    }
};

#endif
