#ifndef FOURCHAMBERCOMMON_H
#define FOURCHAMBERCOMMON_H

#include <QString>

struct FourChamberSubfolders {
    QString SEG, MESH, UVC, UVC_LA, UVC_RA, AFIB, PRESIM, SIMS;
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

/// @brief Meshtool3D static libraries parameters structures
struct M3DParameters {
    QString seg_dir, seg_name, out_dir, out_name;

    bool mesh_from_segmentation, boundary_relabelling, verbose, eval_thickness;
    bool  out_medit, out_carp, out_carp_binary, out_vtk, out_vtk_binary, out_potential;

    double facet_angle, facet_size, facet_distance;
    double cell_rad_edge_ratio, cell_size;
    double rescaleFactor, abs_tol, rel_tol;
    int itr_max;
    int dimKrilovSp;

    M3DParameters()
    :seg_dir(""),
     seg_name(""),
     mesh_from_segmentation(true),
     boundary_relabelling(false),
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
     verbose(true),
     eval_thickness(false),
     out_dir(""),
     out_name(""),
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
    double alpha_endo, alpha_epi, beta_endo, beta_epi;
    QString bad_elem, type, apex_to_base, epi, lv, rv;
    QString directory;
    int nonmyo;

    VFibresParams()
    :alpha_endo(60),
     alpha_epi(-60),
     beta_endo(-65), 
     beta_epi(25), 
     bad_elem(""), 
     type("biv"), 
     apex_to_base(""), 
     epi(""),
     lv(""),
     rv("")
     {}

     inline QString a_endo() { return QString::number(alpha_endo); }; 
     inline QString a_epi() { return QString::number(alpha_epi); }; 
     inline QString b_endo() { return QString::number(beta_endo); }; 
     inline QString b_epi() { return QString::number(beta_epi); };
     inline QString type() { return type; };
     inline QString bad_elem(QString dir="") { return dir + "/" + bad_elem; };
     inline QString apex_to_base(QString dir="") { return dir + "/" + apex_to_base; };
     inline QString epi(QString dir="") { return dir + "/" + epi; };
     inline QString lv(QString dir="") { return dir + "/" + lv; };
     inline QString rv(QString dir="") { return dir + "/" + rv; };
}

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