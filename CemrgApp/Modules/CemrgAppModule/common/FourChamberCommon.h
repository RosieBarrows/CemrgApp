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