#!/usr/bin/env python3

__email__    = "conrad.shyu@nih.gov"
__author__   = "Conrad Shyu"
__version__  = "2.3"
__branch__   = "Bioinformatics and Computational Biosciences Branch"
__company__  = "National Institute of Allergy and Infectious Diseases"
__office__   = "5601 Fishers Lane, Rockville, MD 20852"
__update__   = "3/9/2021"
__project__  = "METAGENOTE"
__modified__ = "Los Alamos National Laboratory, Biosecurity and Public Health"
__contributor__ = "Chienchi Lo (chienchi@lanl.gov)"

import os
import sys
import json
import datetime
import argparse
import subprocess
import xml.dom.minidom
import xml.etree.ElementTree as ET

class SubmitXML():
    """ generate the XML file for submission to NCBI """
    def __init__(self, org = "LANL BPH"):
        """ generate the XML file for SRA submission  """
        self.root = ET.Element('Submission')
        self.library = {}   # library details
        self.fastq = {}     # list of fastq files
        self.sample = {}    # sample attributes
        self.project = {}   # project details
        self.org = org      # center name
        self.package = ""   # package name

    def ReadFile(self, xf):
        """ read the contents of a file """
        lines = []
        with open(xf, "r") as f:
            for r in f.readlines():
                if len(r.strip()) is 0:
                    continue
                lines.append(r.strip())
        return(lines)

    def GetFASTQ(self, fq):
        """ list of FASTQ files and associated attributes; sample name as key, fastq files as list """
        fastq = {}
        for r in self.ReadFile(fq):
            fq = r.split(",")
            fastq[(fq[0]).strip()] = [q.strip() for q in fq[1:]]
        return(fastq)

    def GetAttribute(self, sg):
        """ sample metadata; sample name as key, attributes as dictionary """
        attr = {}
        lines = self.ReadFile(sg)
        keys = [k.strip() for k in (lines[1]).split("\t")]
        for r in lines[2:]:
            v = r.split("\t")
            attr[(v[0]).strip()] = dict(zip(keys[1:], [i.strip() for i in v[1:]]))
        return(attr)
        
    def GetLibraryAttribute(self, exp):
        """ sample library metadata; sample name as key, attributes as dictionary """
        attr = {}
        lines = self.ReadFile(exp)
        keys = [k.strip() for k in (lines[0]).split("\t")]
        for r in lines[1:]:
            v = r.split("\t")
            attr[(v[0]).strip()] = dict(zip(keys[1:], [i.strip() for i in v[1:]]))
        return(attr)
        
    def GetProject(self, pj):
        """ project information """
        proj = {}
        for r in self.ReadFile(pj):
            p = r.split("\t")
            proj[(p[0]).strip()] = (p[1]).strip()
        return(proj)

    def GetPackage(self, sg):
        """ return the package name """
        return((self.ReadFile(sg)[0]).strip("# \n"))

    def SetElement(self, tag, text = "", attr = {}, child = []):
        """ set the basic XML element """
        xt = ET.Element(tag, attrib = attr)
        xt.text = text

        for i in child:
            xt.append(i)

        return(xt)      # end of SetElement

    def SetDescription(self, email, first = "NIAID", last = "BCBB",
        org = "Los Alamos National Laboratory, Biosecurity and Public Health",
        date = datetime.datetime.now().strftime("%Y-%m-%d")):
        """ set the submission description """
        self.org = org      # reset the center name

        self.root.append(self.SetElement("Description", "", {}, [
            self.SetElement("Comment", "New submission. Generated by METAGENOTE on %s" % datetime.datetime.now().strftime("%A %B %d %Y %X")),
            self.SetElement("Organization", "", {"role": "owner", "type": "center"}, [
                self.SetElement("Name", self.org),
                self.SetElement("Contact", "", {"email": email}, [
                    self.SetElement("Name", "", {}, [
                        self.SetElement("First", first),
                        self.SetElement("Last", last)])])
                ]),
            self.SetElement("Hold", "", {"release_date": date})
            ]))

        return(True)    # end of SetDescription

    def SetDescriptor(self, title, desc, url):
        """ auxiliary for external resources """
        xm = []
        xm.append(self.SetElement("Title", title))
        xm.append(self.SetElement("Description", "", {}, [self.SetElement(
            "p", desc + " This submission was powered by METAGENOTE (https://metagenote.niaid.nih.gov).")]))

        for k, v in (json.loads(url)).items():
            xm.append(self.SetElement("ExternalLink", "", {"label": k}, [self.SetElement("URL", v)]))

        return(xm)      # end of SetDescriptor

    def SetBioProject(self, title, name, desc, data, organism, scope,
        url = '{"EDGE COVID19": "https://edge-covid19.edgebioinformatics.org/"}'):
        """ set the bioproject records """
        pts = [self.SetElement("IntendedDataTypeSet", "", {}, [self.SetElement("DataType", data)])]
        if organism:
            pts = pts + [self.SetElement("Organism", "",{},[self.SetElement("OrganismName",organism)])]
        self.root.append(self.SetElement("Action", "", {}, [
            self.SetElement("AddData", "", {"target_db": "BioProject"}, [
                self.SetElement("Data", "", {"content_type": "XML"}, [
                    self.SetElement("XmlContent", "", {}, [
                        self.SetElement("Project", "", {"schema_version": "2.0"}, [
                            self.SetElement("ProjectID", "", {}, [
                                self.SetElement("SPUID", name, {"spuid_namespace": self.org})]),
                            self.SetElement("Descriptor", "", {},
                                self.SetDescriptor(title, desc, url)),      # title, description and external resources
                            self.SetElement("ProjectType", "", {}, [
                                self.SetElement("ProjectTypeSubmission", "", {"sample_scope": scope}, pts)])
                        ])])]),
                self.SetElement("Identifier", "", {}, [
                    self.SetElement("SPUID", name, {"spuid_namespace": self.org})])
                ])]))

        return(True)    # end of SetBioProject

    def SetBioSample(self, title, spuid, sid, name, pkg, attr = {}):
        """ biosample section """
        rf = [
            self.SetElement("SampleId", "", {}, [
                self.SetElement("SPUID", sid , {"spuid_namespace": self.org})]),
            self.SetElement("Descriptor", "", {}, [
                self.SetElement("Title", title),
                self.SetElement("Description", "", {}, [
                    self.SetElement("p", spuid)])
                ]),
            self.SetElement("Organism", "", {}, [
                self.SetElement("OrganismName", name)]),
            self.SetElement("Package", pkg),
            self.SetElement("Attributes", "", {}, [
                self.SetElement("Attribute", attr[k], {"attribute_name": k}) for k in sorted(attr.keys())])]

        if "bioproject_accession" in attr:
            a = [self.SetElement("BioProject", "", {}, [
                self.SetElement("PrimaryId", attr["bioproject_accession"], {"db": "BioProject"})])]
            rf = rf[:3] + a + rf[3:]

        self.root.append(self.SetElement("Action", "", {}, [
            self.SetElement("AddData", "", {"target_db": "BioSample"}, [
                self.SetElement("Data", "", {"content_type": "XML"}, [
                    self.SetElement("XmlContent", "", {}, [
                        self.SetElement("BioSample", "", {"schema_version": "2.0"}, rf)])]),
                self.SetElement("Identifier", "", {}, [
                    self.SetElement("SPUID", sid , {"spuid_namespace": self.org})])
            ])]))

        return(True)    # end of SetBiosample

    def SetFASTQ(self, fastq, title, sid, attr, date, bioprojectID = None):
        """ FASTQ files and sequencing instrument """
        # list of FASTQ files
        fq = [self.SetElement("File", "", {"file_path": f}, [
            self.SetElement("DataType", "generic-data")]) for f in fastq]
            
        if "design_description" in attr and not attr["design_description"]:
            del attr["design_description"]
        if "title" in attr and not attr["title"]:
            del attr["title"]
        for i in ['library_ID','filetype','filename','filename2','filename3','filename4','assembly','fasta_file']:
            attr.pop(i, None)
        
        # instrument and library attributes
        attr["library_name"] = "%s.%s" % (fastq[0], sid)
        at = [self.SetElement("Attribute", attr[k], {"name": k}) for k in attr.keys()]

        if bioprojectID:
            bpid_attr = [self.SetElement("PrimaryId", bioprojectID, {"db": "BioProject"})]
        else:
            bpid_attr = [self.SetElement("SPUID", title, {"spuid_namespace": self.org})]
        # supplemental information
        rf = [
            self.SetElement("AttributeRefId", "", {"name": "BioProject"}, [
                self.SetElement("RefId", "", {}, bpid_attr)]),
            self.SetElement("AttributeRefId", "", {"name": "BioSample"}, [
                self.SetElement("RefId", "", {}, [
                    self.SetElement("SPUID", sid + ":" + date, {"spuid_namespace": self.org})])]),
            self.SetElement("Identifier", "", {}, [
                self.SetElement("LocalId", "%s.%s" % (fastq[0], sid))])]

        self.root.append(self.SetElement("Action", "", {}, [
            self.SetElement("AddFiles", "", {"target_db": "SRA"}, fq + at + rf)]))

        return(True)    # end of SetFASTQ

    def Print(self):
        print((xml.dom.minidom.parseString(ET.tostring(self.root))).toprettyxml(indent = "    ", newl = "\n"))

    def Write(self, name = "submission.xml"):
        #(ET.ElementTree(root)).write("submission.xml", encoding="UTF-8", method="xml")
        with open(name, "w") as to_xml:
            to_xml.write((xml.dom.minidom.parseString(ET.tostring(self.root))).toprettyxml(indent = "    ", newl = "\n"))

    def Run(self, lb, datatype, scope, lf = "listfile.txt", pj = "project.txt", sg = "samplegroup.txt", exp = "sra_experiments.txt"):
        """ main driver function """
        library = lb
        fastq = self.GetFASTQ(lf)           # list of FASTQ files
        sample = self.GetAttribute(sg)      # sample attributes
        project = self.GetProject(pj)       # project details
        package = self.GetPackage(sg)       # package name
        if (os.path.exists(exp)):
            library = self.GetLibraryAttribute(exp)

        k = next(iter(sample))
        library["library_construction_protocol"] = (sample[k])["lib_const_meth"] if "lib_const_meth" in sample[k] else "missing"
        #library["design_description"] = (sample[k])["lib_const_meth"] if "lib_const_meth" in sample[k] else "missing"

        self.SetDescription(
            project["Email"],               # email address
            project["First"],               # first name
            project["Last"],                # last name
            project["Center"],              # center name
            project["Release"])             # hold release date

        if "bioproject_accession" not in sample[k]:
            organism_name = project["Organism"] if "Organism" in project else None
            self.SetBioProject(
                project["ProjectTitle"],    # project title
                project["ProjectName"],     # project name
                project["Description"],     # project description
                datatype,        # project datatype
                organism_name,              # project organism for single and multiisolate projects
                scope,    # sample scope 
                project["Resource"])		# external resources
                
        date = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        # set biosample section
        for k, v in sample.items():
            sid = k + ":" + date
            name = ""
            sample_title= " "
            if 'sample_title' in v:
                if v["sample_title"]:
                    sample_title = v["sample_title"]
                else:
                    del v["sample_title"]
                
            if "organism" in v:
                name = v["organism"]
                del v["organism"]           # not appear as an attribute
            else:
                name = "not applicable"

            self.SetBioSample(sample_title, v["description"], sid, name, package, v)

        # set fastq file section
        fqs = []
        for k, v in fastq.items():
            sid = k
            fqs += [False if f == "null" else True for f in v]
            exp_library = library[k] if os.path.exists(exp) else library
            if "bioproject_accession" not in sample[k]:
                self.SetFASTQ(v, project["ProjectName"], sid, exp_library, date)
            else:
                self.SetFASTQ(v, None, sid, exp_library, date, sample[k]["bioproject_accession"])

        return(not (False in fqs))

def SaveEmail(pf):
    """ auxiliary functions """
    # first, load the file into memory
    with open(pf, "r") as rf:
        lines = rf.readlines()

    # then, write the contents back without the last line
    with open(pf, "w") as wf:
        ue = (lines[len(lines) - 1]).strip()
        wf.writelines([item for item in lines[:-1]])

    # write the email address to a file
    with open(os.path.join(os.path.dirname(pf), ".user-email"), "w") as ef:
        ef.write(ue)

    return(ue)

def CheckStatus(f):
    """ check the status of validate.xml """
    t = ET.ElementTree(file = f)
    s = {"processed-error": False, "failed": False, "processed-ok": True, "submitted": False, "processing": False}
    #return(s["%s" % (t.getroot()).find("Action").attrib["status"]])
    return((t.getroot()).find("Action").attrib["status"])

def ValidateLog(f):
    """ output the validate.xml content to logfile """
    msg = {}
    t = (ET.ElementTree(file = f)).getroot()
    """
    with open(f, "w") as fl:
        fl.write("message: %s\n" % t.find("Action/Response/Message").text)
        fl.write("sample_id: %s\n" % t.find("Action/Response/Object").attrib["spuid"])
    """
    s = {"processed-error": "No", "failed": "No", "processed-ok": "Yes", "submitted": "No", "processing": "No"}
    msg["Passed"] = s[t.find("Action").attrib["status"]]
    if t.find("Action/Response/Message") is not None:
        msg["Report Text"] = t.find("Action/Response/Message").text
    print(json.dumps(msg))
    return(0 if msg["Passed"] == "Yes" else 1)

def main(argv):
    validate_xml = os.path.join(os.path.dirname(argv.project), "validate.xml")
    submit_xml = os.path.join(argv.input_dir, "submission.xml")
    print("Generate Submission.xml file")
    #SaveEmail(args.project) # save the user email address

    ncbi = SubmitXML()      # generate the XML file for SRA submission
    r = ncbi.Run({
        "library_selection": argv.libselection,
        "library_strategy": argv.libstrategy,
        "library_layout": argv.liblayout,
        "library_source": argv.libsource,
        "instrument_platform": argv.platform,
        "instrument_model": argv.instrument,
        "design_description":argv.libdesign,
        "title":argv.libtitle },
        argv.datatype.lower(), argv.samplescope, argv.listfile, argv.project, argv.samplegroup, argv.experiments)

    if (not r):
        print("{'Passed': 'No', 'Report Text': 'Missing sequencing file(s)'}")
        return(1)

    ncbi.Write(submit_xml)
    print("Validate Submission.xml by NCBI ")
    with open(validate_xml, "w") as rp:
        subprocess.call(["curl", "-X", "POST", "-d", "@%s" % submit_xml,
            "https://www.ncbi.nlm.nih.gov/projects/biosample/validate/"], stdout = rp)

    #return(CheckStatus(validate_xml))
    return(ValidateLog(validate_xml))

if __name__ == '__main__':
    argv = argparse.ArgumentParser(description = "Generate XML file for submission to NCBI SRA.",
        formatter_class = argparse.RawTextHelpFormatter)
    argv.add_argument('-f', '--input-dir', dest = 'input_dir', required = True,
        help = 'Diretory where all of your FASTQ files are located.')
    argv.add_argument('-l', '--listfile', type = str, dest = 'listfile', default = "listfile.txt",
        help = 'Name of the list file containing your sample ID to file mapping. Default: listfile.txt.')
    #argv.add_argument('-e', '--email', dest = 'email', default = "metagenote@mail.nih.gov",
    #    help = 'Eamil address for submission notification. Default: metagenote@mail.nih.gov.')
    argv.add_argument('-p', '--projectfile', type = str, dest = 'project', default = "project.txt",
        help = 'Name of the file describing your project. Default: project.txt.')
    argv.add_argument('-m', '--metadatafile', type = str, dest = 'samplegroup', default = "samplegroup.txt",
        help = 'Name of the mimark metadata file. Default: samplegroup.txt.')
    argv.add_argument('-e',  '--experimentfile', dest = 'experiments', default = "sra_experiments.txt", 
        help = 'Name of the sra experments (library and platform). Default: sra_experiments.txt')

    # the following command line parameters were added by Conrad Shyu, July 3, 2017
    argv.add_argument('-a', '--platform', dest = 'platform', default = 'Illumina',
        help = 'Specify platform you are using for sequencing. Default: illumina.')
    #argv.add_argument('-o', '--orientation', dest = 'orientation', default = 'forward',
        #help = 'Specify sequence orientation. Default: forward.')
    argv.add_argument('-s', '--libstrategy', dest = 'libstrategy', default = 'amplicon',
        help = 'Specify library strategy. Default: amplicon.')
    argv.add_argument('-t', '--datatype', dest = 'datatype', default = 'metagenome',
        help = 'Specify datatype. Default: metagenome.')
    argv.add_argument('-r', '--libsource', dest = 'libsource', default = 'metagenomic',
        help = 'Specify library source. Default: metagenomic.')
    argv.add_argument('-n', '--instrument', dest = 'instrument', default = 'Illumina MiSeq',
        help = "Specify instrument. Default: Illumina MiSeq.")
    argv.add_argument('-b', '--libselection', dest = 'libselection', default = 'PCR',
        help = "Specify library selection. Default: PCR.")
    argv.add_argument('-y', '--liblayout', dest = 'liblayout', default = 'paired',
        help = "Specify library layout. Default: paired.")
    argv.add_argument('--libdesign', dest = 'libdesign', default = '',
        help = "Brief description of the methods used to create the sequencing library.")
    argv.add_argument('--libtitle', dest = 'libtitle', default = '',
        help = "Short description that will identify the dataset on public pages. {methodology} of {organism}: {sample info}")
    argv.add_argument('-c', '--samplescope', dest = 'samplescope', default = 'eMultiisolate',
        help = "Specify The scope and purity of the biological sample used for the study. Default: eMultiisolate.")
    

    sys.exit(main(argv.parse_args()))
