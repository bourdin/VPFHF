#include "petsc.h"


#undef __FUNCT__
#define __FUNCT__ "XDMFuniformgridInitialize"
/*
  XDMFuniformgridInitialize: 

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode XDMFuniformgridInitialize(PetscViewer viewer,PetscReal time,const char gridname[])
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"<?xml version=\"1.0\" ?>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t<Domain>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t<Grid Name=\"%s\" GridType=\"Uniform\">\n",gridname);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t<Time Value=\"%f\" />\n",time);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFtopologyAdd"
/*
  XDMFtopologyAdd: 

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode XDMFtopologyAdd(PetscViewer viewer,PetscInt nx,PetscInt ny,PetscInt nz,const char h5filename[],const char coordname[])
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t<Topology TopologyType=\"3DSMesh\" Dimensions=\"%i %i %i\"/>\n",nz,ny,nx);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t<Geometry GeometryType=\"XYZ\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t\t<DataItem Dimensions=\"%i %i %i 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:%s</DataItem>\n",nz,ny,nx,h5filename,coordname);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t</Geometry>\n");CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFattributeAdd"
/*
  XDMFattributeAdd: 

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode XDMFattributeAdd(PetscViewer viewer,PetscInt nx,PetscInt ny,PetscInt nz,PetscInt nfields,const char fieldtype [],const char location[],const char h5filename[],const char fieldname[])
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t<Attribute Name=\"%s\" Active=\"1\" AttributeType=\"%s\" Center=\"%s\">\n",fieldname,fieldtype,location);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t\t<DataItem Dimensions=\"%i %i %i %i\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:%s</DataItem>\n",nz,ny,nx,nfields,h5filename,fieldname);CHKERRQ(ierr);      
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t</Attribute>\n");CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFuniformgridFinalize"
/*
  XDMFuniformgridFinalize: 

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode XDMFuniformgridFinalize(PetscViewer viewer)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t</Grid>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t</Domain>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"</Xdmf>\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "XDMFmultistepInitialize"
/*
  XDMFmultistepInitialize: 

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode XDMFmultistepInitialize(PetscViewer viewer)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"<?xml version=\"1.0\" ?>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t<Domain>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFmultistepAddstep"
/*
  XDMFmultistepAddstep: 

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode XDMFmultistepAddstep(PetscViewer viewer,const char filename[])
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t\t<xi:include href=\"%s\" xpointer=\"xpointer(//Xdmf/Domain/Grid)\" />\n",filename);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFmultistepFinalize"
/*
  XDMFmultistepFinalize: 

  (c) 2010 Blaise Bourdin bourdin@lsu.edu
*/
extern PetscErrorCode XDMFmultistepFinalize(PetscViewer viewer)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"\t\t</Grid>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\t</Domain>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"</Xdmf>\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
