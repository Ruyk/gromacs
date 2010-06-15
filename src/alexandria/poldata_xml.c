/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: poldata_xml.c,v 1.19 2009/05/03 14:19:26 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "gmx_fatal.h"
#include "macros.h"
#include "grompp.h"
#include "smalloc.h"
#include "poldata_xml.h"
#include "xml_util.h"
#include "futil.h"

extern int xmlDoValidityCheckingDefaultValue;
	
#define NN(x) (NULL != (x))

static const char *xmltypes[] = { 
    NULL, 
    "XML_ELEMENT_NODE",
    "XML_ATTRIBUTE_NODE",
    "XML_TEXT_NODE",
    "XML_CDATA_SECTION_NODE",
    "XML_ENTITY_REF_NODE",
    "XML_ENTITY_NODE",
    "XML_PI_NODE",
    "XML_COMMENT_NODE",
    "XML_DOCUMENT_NODE",
    "XML_DOCUMENT_TYPE_NODE",
    "XML_DOCUMENT_FRAG_NODE",
    "XML_NOTATION_NODE",
    "XML_HTML_DOCUMENT_NODE",
    "XML_DTD_NODE",
    "XML_ELEMENT_DECL",
    "XML_ATTRIBUTE_DECL",
    "XML_ENTITY_DECL",
    "XML_NAMESPACE_DECL",
    "XML_XINCLUDE_START",
    "XML_XINCLUDE_END"
};
#define NXMLTYPES asize(xmltypes)
	
enum { 
    exmlGENTOP,
    exmlSMATOMS, exmlPOLAR_UNIT, 
    exmlSMATOM, exmlELEM, exmlNAME, exmlDESC,
    exmlSMNAME, exmlSMTYPE, exmlMILLER_EQUIV, exmlCHARGE,
    exmlNEIGHBORS, 
    exmlGEOMETRY, exmlNUMBONDS, exmlPOLARIZABILITY, exmlSIGPOL, exmlSPREF,
    exmlSMBONDS, exmlLENGTH_UNIT, exmlSMBOND,
    exmlATOM1, exmlATOM2, exmlLENGTH, exmlBONDORDER,
    exmlSMANGLES, exmlANGLE_UNIT, exmlSMANGLE,
    exmlATOM3, exmlANGLE, 
    exmlBSATOMS, exmlBSATOM,
    exmlMILATOMS, exmlTAU_UNIT, exmlAHP_UNIT,
    exmlMILATOM, exmlMILNAME, exmlSPOEL_EQUIV,
    exmlATOMNUMBER, exmlTAU_AHC, exmlALPHA_AHP,
    exmlSYMMETRIC_CHARGES, exmlSYM_CHARGE,
    exmlCENTRAL, exmlATTACHED, exmlNUMATTACH,
    exmlEEMPROPS, exmlEEMPROP, exmlMODEL, exmlJ0, exmlCHI0, exmlZETA, exmlROW,
    exmlEEMPROP_REF, exmlEPREF,
    exmlNR 
};
  
static const char *exml_names[exmlNR] = {
    "gentop",
    "smatoms", "polarizability_unit",
    "smatom", "elem", "name", "description",
    "smname", "smtype", "miller_equiv", "charge",
    "neighbors", 
    "geometry", "numbonds", "polarizability", "sigma_pol", "spref",
    "smbonds", "length_unit", "smbond",
    "atom1", "atom2", "length", "bondorder",
    "smangles", "angle_unit", "smangle",
    "atom3", "angle",
    "bsatoms", "bsatom",
    "milatoms", "tau_ahc_unit", "alpha_ahp_unit",
    "milatom", "milname", "spoel_equiv",
    "atomnumber", "tau_ahc", "alpha_ahp",
    "symmetric_charges", "sym_charge",
    "central", "attached", "numattach",
    "eemprops", "eemprop", "model", "jaa", "chi", "zeta", "row",
    "eemprop_ref", "epref"
};

static void sp(int n, char buf[], int maxindent)
{
    int i;
    if(n>=maxindent)
        n=maxindent-1;
  
    /* Don't indent more than maxindent characters */
    for(i=0; (i<n); i++)
        buf[i] = ' ';
    buf[i] = '\0';
}

static void process_attr(FILE *fp,xmlAttrPtr attr,int elem,
                         int indent,gmx_poldata_t pd,gmx_atomprop_t aps)
{
    char *attrname,*attrval;
    char buf[100];
    int  i,kkk;
    char *xbuf[exmlNR];
  
    for(i=0; (i<exmlNR); i++)
        xbuf[i] = NULL;    
    while (attr != NULL) 
    {
        attrname = (char *)attr->name;
        attrval  = (char *)attr->children->content;
    
#define atest(s) ((strcasecmp(attrname,s) == 0) && (attrval != NULL))
        kkk = find_elem(attrname,exmlNR,exml_names);
        if (-1 != kkk) 
        {
            if (attrval != NULL)
                xbuf[kkk] = strdup(attrval);
            
            if (NULL != fp) {
                sp(indent,buf,99);
                fprintf(fp,"%sProperty: '%s' Value: '%s'\n",buf,attrname,attrval);
            }
        }
        else
        {
            fprintf(stderr,"Ignoring invalid attribute %s\n",attrname);
        }
        attr = attr->next;
#undef atest
    }
    /* Done processing attributes for this element. Let's see if we still need
     *  to interpret them.
     */
     
    switch (elem) {
    case exmlSMATOMS:
        if (NN(xbuf[exmlPOLAR_UNIT])) 
            gmx_poldata_set_spoel_unit(pd,xbuf[exmlPOLAR_UNIT]);
        break;
    case exmlBSATOMS:
        if (NN(xbuf[exmlPOLAR_UNIT])) 
            gmx_poldata_set_bosque_unit(pd,xbuf[exmlPOLAR_UNIT]);
        break;
    case exmlSMANGLES:
        if (NN(xbuf[exmlANGLE_UNIT])) 
            gmx_poldata_set_angle_unit(pd,xbuf[exmlANGLE_UNIT]);
        break;
    case exmlSMBONDS:
        if (NN(xbuf[exmlLENGTH_UNIT])) 
            gmx_poldata_set_length_unit(pd,xbuf[exmlLENGTH_UNIT]);
        break;
    case exmlMILATOMS:
        if (NN(xbuf[exmlTAU_UNIT]) && NN(xbuf[exmlAHP_UNIT]))
            gmx_poldata_set_miller_units(pd,xbuf[exmlTAU_UNIT],xbuf[exmlAHP_UNIT]);
        break;
    case exmlSMATOM:
        if (NN(xbuf[exmlELEM]) && NN(xbuf[exmlMILLER_EQUIV]) && 
            NN(xbuf[exmlCHARGE]) && NN(xbuf[exmlGEOMETRY]) && 
            NN(xbuf[exmlPOLARIZABILITY]) && NN(xbuf[exmlNUMBONDS]) && 
            NN(xbuf[exmlNEIGHBORS]) &&
            NN(xbuf[exmlSPREF]) && NN(xbuf[exmlSMTYPE]) && NN(xbuf[exmlSMNAME]))
            gmx_poldata_add_spoel(pd,xbuf[exmlELEM],
                                  xbuf[exmlDESC] ? xbuf[exmlDESC] : (char *) "",
                                  xbuf[exmlSMNAME],xbuf[exmlSMTYPE],
                                  xbuf[exmlMILLER_EQUIV],
                                  atoi(xbuf[exmlCHARGE]),
                                  xbuf[exmlGEOMETRY],
                                  atoi(xbuf[exmlNUMBONDS]),
                                  xbuf[exmlNEIGHBORS],
                                  atof(xbuf[exmlPOLARIZABILITY]),
                                  NN(xbuf[exmlSIGPOL]) ? atof(xbuf[exmlSIGPOL]) : 0,
                                  xbuf[exmlSPREF]);
        break;
    case exmlMILATOM:
        if (NN(xbuf[exmlMILNAME]) && NN(xbuf[exmlATOMNUMBER]) && 
            NN(xbuf[exmlTAU_AHC]) && NN(xbuf[exmlALPHA_AHP])) 
            gmx_poldata_add_miller(pd,xbuf[exmlMILNAME],atoi(xbuf[exmlATOMNUMBER]),
                                   atof(xbuf[exmlTAU_AHC]),atof(xbuf[exmlALPHA_AHP]),
                                   xbuf[exmlSPOEL_EQUIV]);
        break;
    case exmlBSATOM:
        if (NN(xbuf[exmlELEM]) && NN(xbuf[exmlPOLARIZABILITY]))
            gmx_poldata_add_bosque(pd,xbuf[exmlELEM],atof(xbuf[exmlPOLARIZABILITY]));
        break;
    case exmlSMBOND:
        if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) && 
            NN(xbuf[exmlLENGTH]) && NN(xbuf[exmlBONDORDER]))
            gmx_poldata_add_smbond(pd,xbuf[exmlATOM1],xbuf[exmlATOM2],
                                   atof(xbuf[exmlLENGTH]),atof(xbuf[exmlBONDORDER]));
        break;
    case exmlANGLE:
        if (NN(xbuf[exmlATOM1]) && NN(xbuf[exmlATOM2]) && 
            NN(xbuf[exmlATOM3]) && NN(xbuf[exmlANGLE]))
            gmx_poldata_add_smangle(pd,xbuf[exmlATOM1],xbuf[exmlATOM2],
                                    xbuf[exmlATOM3],atof(xbuf[exmlANGLE]));
        break;
    case exmlSYM_CHARGE:
        if (NN(xbuf[exmlCENTRAL]) && NN(xbuf[exmlATTACHED]) && 
            NN(xbuf[exmlNUMATTACH]))
            gmx_poldata_add_symcharges(pd,xbuf[exmlCENTRAL],
                                       xbuf[exmlATTACHED],
                                       atoi(xbuf[exmlNUMATTACH]));
        break;
    case exmlEEMPROP:
        if (NN(xbuf[exmlMODEL]) && NN(xbuf[exmlNAME]) && 
            NN(xbuf[exmlCHI0])  && NN(xbuf[exmlJ0]) && 
            NN(xbuf[exmlZETA])  && NN(xbuf[exmlCHARGE]) && 
            NN(xbuf[exmlROW])) 
            gmx_poldata_set_eemprops(pd,name2eemtype(xbuf[exmlMODEL]),xbuf[exmlNAME],
                                     atof(xbuf[exmlJ0]),atof(xbuf[exmlCHI0]),
                                     xbuf[exmlZETA],xbuf[exmlCHARGE],xbuf[exmlROW]);
        break;
    case exmlEEMPROP_REF:
        if (NN(xbuf[exmlMODEL]) && NN(xbuf[exmlEPREF]))
            gmx_poldata_set_epref(pd,name2eemtype(xbuf[exmlMODEL]),xbuf[exmlEPREF]);
        break;
    default:
        if (NULL != debug) {
            fprintf(debug,"Unknown combination of attributes:\n");
            for(i=0; (i<exmlNR); i++)
                if (xbuf[i] != NULL)
                    fprintf(debug,"%s = %s\n",exml_names[i],xbuf[i]);
        }
    }
    
    /* Clean up */
    for(i=0; (i<exmlNR); i++)
        if (xbuf[i] != NULL)
            sfree(xbuf[i]);

}

static void process_tree(FILE *fp,xmlNodePtr tree,int parent,int indent,
                         gmx_poldata_t pd,gmx_atomprop_t aps)
{
    int elem;
    char          buf[100];
  
    while (tree != NULL) 
    {
        if (fp) 
        {
            if ((tree->type > 0) && (tree->type < NXMLTYPES))
                fprintf(fp,"Node type %s encountered with name %s\n",
                        xmltypes[tree->type],(char *)tree->name);
            else
                fprintf(fp,"Node type %d encountered\n",tree->type);
        }
    
        if (tree->type == XML_ELEMENT_NODE)
        {
            elem = find_elem((char *)tree->name,exmlNR,exml_names);
            if (fp) 
            {
                sp(indent,buf,99);
                fprintf(fp,"%sElement node name %s\n",buf,(char *)tree->name);
            }
            if (-1 != elem) 
            {
                if (elem != exmlGENTOP)
                    process_attr(fp,tree->properties,elem,indent+2,pd,aps);
                
                if (tree->children)
                    process_tree(fp,tree->children,elem,indent+2,pd,aps);
            }
        }
        tree = tree->next;
    }
}

gmx_poldata_t gmx_poldata_read(const char *fn,gmx_atomprop_t aps)
{
    xmlDocPtr     doc;
    int           i,npd;
    gmx_poldata_t pd;
    char *fn2;
  
    if (fn)
        fn2 = (char *)gmxlibfn(fn);
    else
        fn2 = (char *)gmxlibfn("alexandria.ff/gentop.dat");
  
    if (NULL != debug)
        fprintf(debug,"Opening library file %s\n", fn2);

    xmlDoValidityCheckingDefaultValue = 0;
    if ((doc = xmlParseFile(fn2)) == NULL) {
        fprintf(stderr,"\nError reading XML file %s. Run a syntax checker such as nsgmls.\n",
                fn2);
        return NULL;
    }
    pd = gmx_poldata_init();
    process_tree(debug,doc->children,0,0,pd,aps);

    xmlFreeDoc(doc);
  
    if (NULL != debug)
        gmx_poldata_write("pdout.dat",pd,aps,0);
  
    return pd;
}

static void add_xml_poldata(xmlNodePtr parent,gmx_poldata_t pd,
                            gmx_atomprop_t aps)
{
    xmlNodePtr child,grandchild,comp;
    int    i,charge,atomnumber,numbonds,
        numattach,element,model;
    char *elem,*miller_equiv,*geometry,*name,*smtype,*spoel_equiv,*spref,*blu,
        *atom1,*atom2,*atom3,*tmp,*central,*attached,*tau_unit,*ahp_unit,
        *epref,*desc;
    char *neighbors,*zeta,*qstr,*rowstr;
    double polarizability,sig_pol,length,tau_ahc,alpha_ahp,angle,J0,chi0,
        bondorder;
  
    child = add_xml_child(parent,exml_names[exmlSMATOMS]);
    tmp = gmx_poldata_get_spoel_unit(pd);
    if (tmp)
        add_xml_char(child,exml_names[exmlPOLAR_UNIT],tmp);
  
    while ((name = gmx_poldata_get_spoel(pd,NULL,&elem,&desc,&smtype,&miller_equiv,
                                         &charge,&geometry,
                                         &numbonds,&neighbors,
                                         &polarizability,&sig_pol,&spref)) != NULL) {
        grandchild = add_xml_child(child,exml_names[exmlSMATOM]);
        add_xml_char(grandchild,exml_names[exmlELEM],elem);
        add_xml_char(grandchild,exml_names[exmlDESC],desc);
        add_xml_char(grandchild,exml_names[exmlSMNAME],name);
        add_xml_char(grandchild,exml_names[exmlSMTYPE],smtype);
        add_xml_char(grandchild,exml_names[exmlMILLER_EQUIV],miller_equiv);
        add_xml_char(grandchild,exml_names[exmlGEOMETRY],geometry);
        add_xml_int(grandchild,exml_names[exmlNUMBONDS],numbonds);
        add_xml_char(grandchild,exml_names[exmlNEIGHBORS],neighbors);
        add_xml_int(grandchild,exml_names[exmlCHARGE],charge);
        add_xml_double(grandchild,exml_names[exmlPOLARIZABILITY],polarizability);
        add_xml_double(grandchild,exml_names[exmlSIGPOL],sig_pol);
        add_xml_char(grandchild,exml_names[exmlSPREF],spref);
    }

    child = add_xml_child(parent,exml_names[exmlSMBONDS]);
    if ((blu = gmx_poldata_get_length_unit(pd)) != NULL)
        add_xml_char(child,exml_names[exmlLENGTH_UNIT],blu);
    while (gmx_poldata_get_smbond(pd,&atom1,&atom2,&length,&bondorder) == 1) {
        grandchild = add_xml_child(child,exml_names[exmlSMBOND]);
        add_xml_char(grandchild,exml_names[exmlATOM1],atom1);
        add_xml_char(grandchild,exml_names[exmlATOM2],atom2);
        add_xml_double(grandchild,exml_names[exmlLENGTH],length);
        add_xml_double(grandchild,exml_names[exmlBONDORDER],bondorder);
    }
  
    child = add_xml_child(parent,exml_names[exmlSMANGLES]);
    if ((blu = gmx_poldata_get_angle_unit(pd)) != NULL)
        add_xml_char(child,exml_names[exmlANGLE_UNIT],blu);
    while (gmx_poldata_get_smangle(pd,&atom1,&atom2,&atom3,&angle) == 1) {
        grandchild = add_xml_child(child,exml_names[exmlSMANGLE]);
        add_xml_char(grandchild,exml_names[exmlATOM1],atom1);
        add_xml_char(grandchild,exml_names[exmlATOM2],atom2);
        add_xml_char(grandchild,exml_names[exmlATOM3],atom3);
        add_xml_double(grandchild,exml_names[exmlANGLE],angle);
    }
  
    child = add_xml_child(parent,exml_names[exmlBSATOMS]);
    if ((tmp = gmx_poldata_get_bosque_unit(pd)) != NULL)
        add_xml_char(child,exml_names[exmlPOLAR_UNIT],tmp);
  
    while ((name = gmx_poldata_get_bosque(pd,NULL,NULL,&polarizability)) != NULL) {
        grandchild = add_xml_child(child,exml_names[exmlBSATOM]);
        add_xml_char(grandchild,exml_names[exmlELEM],name);
        add_xml_double(grandchild,exml_names[exmlPOLARIZABILITY],polarizability);
    }
    child = add_xml_child(parent,exml_names[exmlMILATOMS]);
    gmx_poldata_get_miller_units(pd,&tau_unit,&ahp_unit);
    if (tau_unit)
        add_xml_char(child,exml_names[exmlTAU_UNIT],tau_unit);
    if (ahp_unit)
        add_xml_char(child,exml_names[exmlAHP_UNIT],ahp_unit);
  
    while ((name = gmx_poldata_get_miller(pd,NULL,&atomnumber,&tau_ahc,&alpha_ahp,&spoel_equiv)) != NULL) {
        grandchild = add_xml_child(child,exml_names[exmlMILATOM]);
        add_xml_char(grandchild,exml_names[exmlMILNAME],name);
        add_xml_int(grandchild,exml_names[exmlATOMNUMBER],atomnumber);
        add_xml_double(grandchild,exml_names[exmlTAU_AHC],tau_ahc);
        add_xml_double(grandchild,exml_names[exmlALPHA_AHP],alpha_ahp);
        if (spoel_equiv)
            add_xml_char(grandchild,exml_names[exmlSPOEL_EQUIV],spoel_equiv);
    }

    child = add_xml_child(parent,exml_names[exmlSYMMETRIC_CHARGES]);
  
    while (gmx_poldata_get_symcharges(pd,&central,&attached,&numattach) == 1) {
        grandchild = add_xml_child(child,exml_names[exmlSYM_CHARGE]);
        add_xml_char(grandchild,exml_names[exmlCENTRAL],central);
        add_xml_char(grandchild,exml_names[exmlATTACHED],attached);
        add_xml_int(grandchild,exml_names[exmlNUMATTACH],numattach);
    }
  
    child = add_xml_child(parent,exml_names[exmlEEMPROPS]);
    while (gmx_poldata_get_eemprops(pd,&model,&name,&J0,&chi0,&zeta,&qstr,&rowstr) == 1) {
        grandchild = add_xml_child(child,exml_names[exmlEEMPROP]);
        add_xml_char(grandchild,exml_names[exmlMODEL],get_eemtype_name(model));
        add_xml_char(grandchild,exml_names[exmlNAME],name);
        add_xml_double(grandchild,exml_names[exmlJ0],J0);
        add_xml_double(grandchild,exml_names[exmlCHI0],chi0);
        add_xml_char(grandchild,exml_names[exmlZETA],zeta);
        add_xml_char(grandchild,exml_names[exmlCHARGE],qstr);
        add_xml_char(grandchild,exml_names[exmlROW],rowstr);
    }
    while (gmx_poldata_list_epref(pd,&model,&epref) == 1) {
        grandchild = add_xml_child(child,exml_names[exmlEEMPROP_REF]);
        add_xml_char(grandchild,exml_names[exmlMODEL],get_eemtype_name(model));
        add_xml_char(grandchild,exml_names[exmlEPREF],epref);
    }
}

void gmx_poldata_write(const char *fn,gmx_poldata_t pd,gmx_atomprop_t aps,
                       int compress)
{
    xmlDocPtr  doc;
    xmlDtdPtr  dtd;
    xmlNodePtr myroot;
    int        i,nmt;
    xmlChar    *libdtdname,*dtdname,*gmx;
  
    gmx        = (xmlChar *) "gentop";
    dtdname    = (xmlChar *) "gentop.dtd";
    libdtdname = dtdname;
  
    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == NULL)
        gmx_fatal(FARGS,"Creating XML document","");
    
    if ((dtd = xmlCreateIntSubset(doc,dtdname,libdtdname,dtdname)) == NULL)
        gmx_fatal(FARGS,"Creating XML DTD","");
    
    if ((myroot = xmlNewDocNode(doc,NULL,gmx,NULL)) == NULL)
        gmx_fatal(FARGS,"Creating root element","");
    dtd->next = myroot;
    myroot->prev = (xmlNodePtr) dtd;
    
    /* Add molecule definitions */
    add_xml_poldata(myroot,pd,aps);

    xmlSetDocCompressMode(doc,compress);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(fn,doc,"ISO-8859-1",2) == 0)
        gmx_fatal(FARGS,"Saving file",fn);
    xmlFreeDoc(doc);
}


