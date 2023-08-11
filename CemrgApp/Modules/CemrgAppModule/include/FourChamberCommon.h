#ifndef FOURCHAMBERCOMMON_H
#define FOURCHAMBERCOMMON_H

#include <QString>
#include <QJsonObject>
#include <QJsonArray>

enum ManualPointsType { CYLINDERS, SLICERS, VALVE_PLAINS };
enum SegmentationTagsType {  
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

struct SegmentationTagsStruct {
    SegmentationTagsType stt;
    SegmentationTagsStruct() : stt(SegmentationTagsType::BACKGROUND) {}
    SegmentationTagsStruct(SegmentationTagsType usrStt) : stt(usrStt) {}

    int GetTag() { return static_cast<int>(stt); } 
    void SetTag(SegmentationTagsType usrStt) { stt = usrStt; }

    void SetTagFromString(QString tagName) { 
        if (tagName.compare("BACKGROUND", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::BACKGROUND); return; }
        if (tagName.compare("BLOODPOOL", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::BLOODPOOL); return; }
        if (tagName.compare("LEFT_VENTRICLE", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::LEFT_VENTRICLE); return; }
        if (tagName.compare("RIGHT_VENTRICLE", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::RIGHT_VENTRICLE); return; }
        if (tagName.compare("LEFT_ATRIUM", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::LEFT_ATRIUM); return; }
        if (tagName.compare("RIGHT_ATRIUM", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::RIGHT_ATRIUM); return; }
        if (tagName.compare("AORTA", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::AORTA); return; }
        if (tagName.compare("PULMONARY_ARTERY", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::PULMONARY_ARTERY); return; }
        if (tagName.compare("LSPV", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::LSPV); return; }
        if (tagName.compare("LIPV", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::LIPV); return; }
        if (tagName.compare("RSPV", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::RSPV); return; }
        if (tagName.compare("RIPV", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::RIPV); return; }
        if (tagName.compare("LAA", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::LAA); return; }
        if (tagName.compare("SVC", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::SVC); return; }
        if (tagName.compare("IVC", Qt::CaseInsensitive) == 0) { SetTag(SegmentationTagsType::IVC); return; }
    }

    QString TagName() {
        switch (stt) {
        case SegmentationTagsType::BACKGROUND: return "BACKGROUND";
        case SegmentationTagsType::BLOODPOOL: return "BLOODPOOL";
        case SegmentationTagsType::LEFT_VENTRICLE: return "LEFT_VENTRICLE";
        case SegmentationTagsType::RIGHT_VENTRICLE: return "RIGHT_VENTRICLE";
        case SegmentationTagsType::LEFT_ATRIUM: return "LEFT_ATRIUM";
        case SegmentationTagsType::RIGHT_ATRIUM: return "RIGHT_ATRIUM";
        case SegmentationTagsType::AORTA: return "AORTA";
        case SegmentationTagsType::PULMONARY_ARTERY: return "PULMONARY_ARTERY";
        case SegmentationTagsType::LSPV: return "LSPV";
        case SegmentationTagsType::LIPV: return "LIPV";
        case SegmentationTagsType::RSPV: return "RSPV";
        case SegmentationTagsType::RIPV: return "RIPV";
        case SegmentationTagsType::LAA: return "LAA";
        case SegmentationTagsType::SVC: return "SVC";
        case SegmentationTagsType::IVC: return "IVC";
        default: return "BACKGROUND";
        }
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

struct SegmentationPointsIds
{
    QStringList cylinders, slicers, valve_plains;
    SegmentationPointsIds(){
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

    QStringList GetPointLabelOptions(ManualPointsType mpt){
        QStringList res = QStringList();

        switch (mpt) {
            case ManualPointsType::CYLINDERS :
                res = CYLINDERS();
                break;

            case ManualPointsType::SLICERS : 
                res = SLICERS();
                break;

            case ManualPointsType::VALVE_PLAINS :
                res = VALVE_PLAINS();
                break;
            default:
                break;
            return res;
        }
        return res;
    }

    QString title(ManualPointsType mpt) {
        QString res = QString();
        switch (mpt) {
            case ManualPointsType::CYLINDERS :
                res = "Cylinders";
                break;

            case ManualPointsType::SLICERS : 
                res = "Slicers";
                break;

            case ManualPointsType::VALVE_PLAINS :
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
    QString _directory, _meshname;
    int _nonmyo;

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
     _rv("")
     {}

     inline QString a_endo() { return QString::number(_alpha_endo); }; 
     inline QString a_epi() { return QString::number(_alpha_epi); }; 
     inline QString b_endo() { return QString::number(_beta_endo); }; 
     inline QString b_epi() { return QString::number(_beta_epi); };
     inline QString type() { return _type; };
     inline QString bad_elem(QString dir="") { return dir + "/" + _bad_elem; };
     inline QString apex_to_base(QString dir="") { return dir + "/" + _apex_to_base; };
     inline QString meshname(QString dir="") { return dir + "/" + _meshname; };
     inline QString epi(QString dir="") { return dir + "/" + _epi; };
     inline QString lv(QString dir="") { return dir + "/" + _lv; };
     inline QString rv(QString dir="") { return dir + "/" + _rv; };
};

/*
arguments << "-m" << directory + "/" + m;
-m, --meshname=STRING      basename of model files

arguments << "--type" << type;

arguments << "-a" << directory + "/" + a;
-a, --apex_to_base=STRING  Solution with 1 at the base an 0 at the apex

arguments << "-e" << directory + "/" + e;
-e, --epi=STRING           Solution with 1 at the epi, 0 at lv/rv endo

arguments << "-l" << directory + "/" + l;
-l, --lv=STRING            Solution with 1 at lv, 0 at epi and rv, needed for

arguments << "-r" << directory + "/" + r;
-r, --rv=STRING            Solution with 1 at rv, 0 at epi and lv, needed for


-o, --output=STRING        filename to write fibres to
-b, --badelem=STRING       filename to write bad element indices to
-t, --type=STRING          type of mesh you are generating fibers for                             (possible values="biv", "lv" default=`biv')
-n, --nonmyo=INT           nonmyocardial tags to ignore


--alpha_endo=FLOAT     Fiber rotation angle on the endocardial surfaces.
--alpha_epi=FLOAT      Fiber rotation angle on the epicardial surfaces.
--beta_endo=FLOAT      Sheet rotation angle on the endocardial surfaces
--beta_epi=FLOAT       Sheet rotation angle on the epicardial surfaces

  -h, --help                 Print help and exit
  -V, --version              Print version and exit
  -m, --meshname=STRING      basename of model files  (default=`')
  -o, --output=STRING        filename to write fibres to (default=`fibres.lon')
  -b, --badelem=STRING       filename to write bad element indices to (default=`')
  -t, --type=STRING          type of mesh you are generating fibers for (possible values="biv", "lv" default=`biv')
  -n, --nonmyo=INT           nonmyocardial tags to ignore

Solutions to Laplace's equation for fiber computation:

  -a, --apex_to_base=STRING  Solution with 1 at the base an 0 at the apex (default=`')
  -e, --epi=STRING           Solution with 1 at the epi, 0 at lv/rv endo (default=`')
  -l, --lv=STRING            Solution with 1 at lv, 0 at epi and rv, needed for biv meshes  (default=`')
  -r, --rv=STRING            Solution with 1 at rv, 0 at epi and lv, needed for biv meshes  (default=`')

Endo/epicardial fiber angles (degrees):

      --alpha_endo=FLOAT     Fiber rotation angle on the endocardial surfaces. (default=`40')
      --alpha_epi=FLOAT      Fiber rotation angle on the epicardial surfaces. (default=`-50')
      --beta_endo=FLOAT      Sheet rotation angle on the endocardial surfaces (default=`-65')
      --beta_epi=FLOAT       Sheet rotation angle on the epicardial surfaces (default=`25')
*/

enum PointsNamesType {
    SVC1, SVC2, SVC3, IVC1, IVC2, IVC3, Ao1, Ao2, Ao3, PArt1, PArt2, PArt3,
    SVC_SLICER1, SVC_SLICER2, SVC_SLICER3, SVC_TIP, IVC_SLICER1, IVC_SLICER2, IVC_SLICER3, IVC_TIP, Ao_TIP, PArt_TIP,
    Ao_WT_TIP, PArt_WT_TIP
};

class BasePointsType {

public:
    BasePointsType(int numEnumValues)
        : points(3 * numEnumValues), pointsSet(false) {}
    
    BasePointsType(ManualPointsType mpt)
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

    inline void PointSet(bool value) { pointsSet = value; };
    inline void PointSetOn() { PointSet(true); };
    inline void PointSetOff() { PointSet(false); };
    inline bool IsPointSet() { return pointsSet; };

protected:
    std::vector<double> points;
    bool pointsSet;

private: 
    int GetNumEnumValues(ManualPointsType mpt) {
        switch (mpt) {
            case ManualPointsType::CYLINDERS:
                return static_cast<int>(PArt3) - static_cast<int>(SVC1) + 1; // Or any calculation that fits your needs
            case ManualPointsType::SLICERS:
                return static_cast<int>(PArt_TIP) - static_cast<int>(SVC_SLICER1) + 1; // Calculate for SlicersPointsType
            case ManualPointsType::VALVE_PLAINS:
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
