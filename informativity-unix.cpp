#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

#include <unistd.h>

using namespace std;

/** struct geneInfo
/// stores information about the blast matches for genome of interest
**/
struct geneInfo{
    string geneID;
    float S1;   //sequence identity against G
    float Q1;   //query coverage against G
    float S2;   //sequence identity against Z
    float Q2;   //query coverage against Z
    float SH;   //sequence identity against S
    float QH;   //query coverage against S
};

/** get_path
/// determines the path working in and the name of the run where all data will be written
**/
bool get_path(string &path, string &run_name)
{
    /** Get folder structure set up **/
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) == NULL)
        perror("getcwd() error");

    path=string(cwd)+"/";

    cout << "Name of run: ";
    string command;
    getline(cin,run_name);
    command="mkdir "+run_name;
    system(command.c_str());
    cout << "Directory " << run_name << " has been created. All output and temporary files will be written to this folder.\n";
    path+=run_name+"/";
    return true;
}

/** file_exists
/// determines if the file indicated by the user/ written by the program exists
**/
inline bool file_exists(const std::string& name)
{
    ifstream f(name.c_str());
    return f.good();
}

/** codon_to_amino_acid
/// return aa 1-letter code for codon
/// returns x for aa if invalid codon or stop codon
**/
char codon_to_amino_acid(char a, char b, char c)
{
	if(a=='t' && b=='t')
	{
		if( (c=='t') || (c=='c') ) return 'f';
		else return 'l';
	}

	if(a=='c' && b=='t') return 'l';

	if(a=='a' && b=='t')
	{
		if(c=='g') return 'm';
		else return 'i';
	}
	if(a=='g' && b=='t') return 'v';

	if(a=='t' && b=='c') return 's';

	if(a=='c' && b=='c') return 'p';

	if(a=='a' && b=='c') return 't';

	if(a=='g' && b=='c') return 'a';

	if(a=='t' && b=='a')
	{
		if( (c=='t') || (c=='c') ) return 'y';
		else return 'x';
	}
	if(a=='c' && b=='a')
	{
		if( (c=='t') || (c=='c') ) return 'h';
		else return 'q';
	}
	if(a=='a' && b=='a')
	{
		if( (c=='t') || (c=='c') ) return 'n';
		return 'k';
	}
	if(a=='g' && b=='a')
	{
		if( (c=='t') || (c=='c') ) return 'd';
		else return 'e';
	}

	if(a=='t' && b=='g')
	{
		if( (c=='t') || (c=='c') ) return 'c';
		if(c=='a') return 'x';
		if(c=='g') return 'w';
	}

	if(a=='c' && b=='g') return 'r';

	if(a=='a' && b=='g')
	{
		if( (c=='t') || (c=='c') ) return 's';
		else return 'r';
	}

	if(a=='g' && b=='g') return 'g';

	return 'x';
}

/** translate
/// translate sequence nt to aa
**/
string translate(string nt_seq)
{
    string aa;
    for(int i=0;i<nt_seq.length();i++)
    {
        if((i+2)<nt_seq.length())
            aa+=codon_to_amino_acid(nt_seq[i],nt_seq[i+1],nt_seq[i+2]);
        i+=2;
    }
    return aa;
}

/** reverse_complement
/// return reverse complement of sequence
**/
string reverse_complement(string nt_seq)
{
    string rc="";
    for(int i=nt_seq.length()-1;i>=0;i--)
    {
        switch(nt_seq[i]){
        case 'a': rc+='t'; break;
        case 't': rc+='a'; break;
        case 'c': rc+='g'; break;
        case 'g': rc+='c'; break;
        }
    }
    return rc;
}

/** generate_s_fasta_for_db
/// takes contigs and generates protein database of all 6 orfs
**/
int generate_s_fasta_for_db(string path, string file)
{
    ifstream in;
    ofstream out;
    string line,header;
    string contig="";
    in.open(file.c_str());
    line=path+file+".6orfs";
    out.open(line.c_str());
    bool first=true;
    while(in.peek()!=EOF)
    {
        getline(in,line,'\n');
        if(line[0]=='>')
        {
            if(first) {out << line << "_ORF1" << endl; header=line; first=false;}
            else
            {
                //calculate
                for(int i=0;i<contig.length();i++) {if(contig[i]>64) contig[i]+=32;}
                out << translate(contig) << endl;  //reading frame 1
                out << header << "_ORF1rc" << endl;
                out << translate(reverse_complement(contig)) << endl;   //reading frame 1'
                contig.erase(contig.begin());
                out << header << "_ORF2" << endl;
                out << translate(contig) << endl;   //reading frame 2
                out << header << "_ORF2rc" << endl;
                out << translate(reverse_complement(contig)) << endl;   //reading frame 2'
                contig.erase(contig.begin());
                out << header << "_ORF3" << endl;
                out << translate(contig) << endl;   //reading frame 3
                out << header << "_ORF3rc" << endl;
                out << translate(reverse_complement(contig)) << endl;   //reading frame 3'
                out << line << "_ORF1" << endl;
                header=line;
                contig="";
            }
        }
        else
            contig+=line;
     }
    for(int i=0;i<contig.length();i++) {if(contig[i]>64) contig[i]+=32;}
    out << translate(contig) << endl;  //reading frame 1
    out << header << "_ORF1rc" << endl;
    out << translate(reverse_complement(contig)) << endl;   //reading frame 1'
    contig.erase(contig.begin());
    out << header << "_ORF2" << endl;
    out << translate(contig) << endl;   //reading frame 2
    out << header << "_ORF2rc" << endl;
    out << translate(reverse_complement(contig)) << endl;   //reading frame 2'
    contig.erase(contig.begin());
    out << header << "_ORF3" << endl;
    out << translate(contig) << endl;   //reading frame 3
    out << header << "_ORF3rc" << endl;
    out << translate(reverse_complement(contig)) << endl;   //reading frame 3'

    in.clear(); in.close();
    out.clear(); out.close();
}

/** create_database
/// creates the three databases
**/
bool create_database(string path, char * file_name, string db_name, string db_type)
{
    if(!file_exists(file_name))
    {
        cout << file_name << " does not exist. Check the path.\nProgram terminating.\n";
        return false;
    }
    string command;
    command="makeblastdb -in " + string(file_name) + " -out " + path + db_name + " -title " + db_name + " -dbtype " + db_type;
    system(command.c_str());
    return true;
}

/** blast
/// performs blast using blast++
**/
bool blast(string blastCommand, string path, string query_file, string db_name, string output_file_name)
{
    if(!file_exists(query_file))
    {
        cout << query_file << " does not exist. Check the path.\nProgram terminating.\n";
        return false;
    }
    string command;
    //output format: qseqid, sseqid, qcovs, pident, length, evalue, bitscore
    command=blastCommand + " -db " + path + db_name + " -evalue 0.0001 -max_hsps 1 -max_target_seqs 1 -query " + query_file + " -outfmt \"10 qseqid sseqid qcovs pident length evalue bitscore\" -out " + path + output_file_name;
    system(command.c_str());
    return true;
}

/** write_matrix
/// reads the blast output files and writes the results to csv
**/
bool write_matrix(string path, geneInfo * result_matrix, int matrix_size)
{
    ifstream in;
    int i,index;
    char comma;
    float score;
    string qseqid,sseqid;
    float qcovs,pident,length,evalue,bitscore;

    string line=path+"xBg.blastx";
    in.open(line.c_str());
    while(in.peek()!=EOF)
    {
        //output format: qseqid, sseqid, qcovs, pident, length, evalue, bitscore
        getline(in,qseqid,','); getline(in,sseqid,','); in>>qcovs; in>>comma; in>>pident; in>>comma;
        in>>length; in>>comma; in>>evalue; in>>comma; in>>bitscore; getline(in,line);
        score=qcovs*pident;
        index=-1;
        for(i=0;i<matrix_size;i++)
        {
            if(result_matrix[i].geneID==qseqid) {index=i; i=matrix_size;}
        }
        if(index>-1)
        {
            if(score>(result_matrix[index].S1*result_matrix[index].Q1)) {result_matrix[index].S1=pident; result_matrix[index].Q1=qcovs;}
        }
    }
    in.clear();
    in.close();

    line=path+"xBz.blastx";
    in.open(line.c_str());
    while(in.peek()!=EOF)
    {
        //output format: qseqid, sseqid, qcovs, pident, length, evalue, bitscore
        getline(in,qseqid,','); getline(in,sseqid,','); in>>qcovs; in>>comma; in>>pident; in>>comma;
        in>>length; in>>comma; in>>evalue; in>>comma; in>>bitscore; getline(in,line);
        score=qcovs*pident;
        index=-1;
        for(i=0;i<matrix_size;i++)
        {
            if(result_matrix[i].geneID==qseqid) {index=i; i=matrix_size;}
        }
        if(index>-1)
        {
            if(score>(result_matrix[index].S2*result_matrix[index].Q2)) {result_matrix[index].S2=pident; result_matrix[index].Q2=qcovs;}
        }
    }
    in.clear();
    in.close();

    line=path+"xBs.blastx";
    if(file_exists(line))
    {
        in.open(line.c_str());
        while(in.peek()!=EOF)
        {
            //output format: qseqid, sseqid, qcovs, pident, length, evalue, bitscore
            getline(in,qseqid,','); getline(in,sseqid,','); in>>qcovs; in>>comma; in>>pident; in>>comma;
            in>>length; in>>comma; in>>evalue; in>>comma; in>>bitscore; getline(in,line);
            score=qcovs*pident;
            index=-1;
            for(i=0;i<matrix_size;i++)
            {
                if(result_matrix[i].geneID==qseqid) {index=i; i=matrix_size;}
            }
            if(index>-1)
            {
                if(score>(result_matrix[index].SH*result_matrix[index].QH)) {result_matrix[index].SH=pident; result_matrix[index].QH=qcovs;}
            }
        }
        in.clear();
        in.close();
    }

    ofstream out;
    line=path+"Informativity_Matrix.csv";
    out.open(line.c_str());
    out << "GeneID,S1,Q1,S2,Q2,SH,QH" << endl;
    for(i=0;i<matrix_size;i++)
    {
        out << result_matrix[i].geneID << "," << result_matrix[i].S1 << "," << result_matrix[i].Q1 << ",";
        out << result_matrix[i].S2 << "," << result_matrix[i].Q2 << ",";
        out << result_matrix[i].SH << "," << result_matrix[i].QH << endl;
    }
    out.clear();
    out.close();
    cout << "Informativity_Matrix.csv has been written to file." << endl;
    return true;
}

/** get_geneIds
/// determines the geneIds generating any hits (against g,z, or s)
**/
vector<string> get_geneIds(string path, int argc)
{
    vector<string> xGenes;
    vector<string>::iterator it;
    ifstream in;
    string line=path+"xBg.blastx";
    in.open(line.c_str());
    while(in.peek()!=EOF)
    {
        getline(in,line,',');
        if(it==xGenes.end()) xGenes.push_back(line);
        getline(in,line);
    }
    in.clear();
    in.close();

    line=path+"xBz.blastx";
    in.open(line.c_str());
    while(in.peek()!=EOF)
    {
        getline(in,line,',');
        it=std::find(xGenes.begin(),xGenes.end(),line);
        if(it==xGenes.end()) xGenes.push_back(line);
        getline(in,line);
    }
    in.clear();
    in.close();

    if(argc==6)
    {
        line=path+"xBs.blastx";
        in.open(line.c_str());
        while(in.peek()!=EOF)
        {
            getline(in,line,',');
            it=std::find(xGenes.begin(),xGenes.end(),line);
            if(it==xGenes.end()) xGenes.push_back(line);
            getline(in,line);
        }
        in.clear();
        in.close();
    }
    std::sort(xGenes.begin(),xGenes.end());
    return xGenes;
}


int main(int argc, char* argv[])
{
    /*** [1]=xNT file; [2]=xAA; [3]=gAA; [4]=zAA; [5]=sAA #5 is option if no sample argc=5***/
    if(argc<5) {cout << "Parameters not entered correctly\nRefer to README.\n"; return 0;}

    string path,run_name;
    if(!get_path(path,run_name)) return 0;

    //create databases
    if(!create_database(path,argv[3],"gDB","prot")) return 0;    //make g database
    if(!create_database(path,argv[4],"zDB","prot")) return 0;    //make z database
    if(argc==6)
    {
        generate_s_fasta_for_db(path,argv[5]);
        string line=path+string(argv[5])+".6orfs";
        if(!create_database(path,const_cast<char*>(line.c_str()),"sDB","prot"))
            return 0;
    }  //make s database

    if(!blast("blastx",path,argv[1],"gDB","xBg.blastx")) return 0;  //blast x against g
    if(!blast("blastx",path,argv[1],"zDB","xBz.blastx")) return 0;  //blast x against z
    if(argc==6) {if(!blast("blastx",path,argv[1],"sDB","xBs.blastx")) return 0;}    //blast x against s

    vector<string> geneIDs=get_geneIds(path, argc); //get geneIDs
    geneInfo * result_matrix=new geneInfo[geneIDs.size()];
    for(int i=0;i<geneIDs.size();i++)
    {
        result_matrix[i].geneID=geneIDs[i];
        result_matrix[i].S1=result_matrix[i].Q1=result_matrix[i].S2=result_matrix[i].Q2=result_matrix[i].SH=result_matrix[i].QH=0;
    }
    write_matrix(path,result_matrix,geneIDs.size());
    geneIDs.clear();

    return 0;
}
