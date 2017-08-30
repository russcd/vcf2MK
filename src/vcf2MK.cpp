#include <string.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <regex>
#include <ctype.h>
#include <time.h>
#include "codon_counts.h"
#include "fold_counts.h" 

using namespace std;

/*

 This program will compute MK stats from a given vcf file and outgroup which can be either maf of vcf but must be aligned to the same reference sequence as the focal population vcf file
 
	compiled on odyssey cluster (gcc 4.7.7) :

			g++ -o vcf2MK -std=c++0x vcf2MK.cpp
 
*/

class cmd_line {
	
public:
	
    /// EXTREMELY BASIC FILTERING FOR VCF AND MAF FILES
    float minimum_focal_vcf_quality ;
    float minimum_outgroup_quality ;           /// either score if maf or qual if vcf

    bool output_invariant ;
    
    bool accept_lower ;
    bool annotate_fold ;
    
    string outgroup_file ;
    string focal_population_vcf_file ;
    string cds_file ;
    
    //// BASIC (VERBOSE) OUTPUT WILL GO TO COUT
    string output_file ;    /// THIS WILL HOLD AGGREGATED MK INFORMATION
    
    void read_cmd_line ( int argc, char *argv[] ) ;
	
} ;

void cmd_line::read_cmd_line ( int argc, char *argv[] ) {
	
    minimum_focal_vcf_quality = 0 ;
    minimum_outgroup_quality = 0 ;
    
    accept_lower = false ;
    annotate_fold = false ;
    output_invariant = true ;

    focal_population_vcf_file = "null" ;
    outgroup_file = "null" ;
    cds_file = "null" ;
    output_file = "MK_counts.txt" ;
    
	for (int i=1; i<argc; i++) {

        if ( strcmp(argv[i],"-f") == 0 ) {
            annotate_fold = true ;
        }
        if ( strcmp(argv[i],"--al") == 0 ) {
            accept_lower = true ; 
        }
        if ( strcmp(argv[i],"-q") == 0 ) {
            minimum_focal_vcf_quality = atof( argv[++i] ) ;
        }
        if ( strcmp(argv[i],"--oq") == 0 ) {
            minimum_outgroup_quality = atof( argv[++i] ) ;
        }
        if ( strcmp(argv[i],"-o") == 0 ) {
            outgroup_file = argv[++i] ;
        }
        if ( strcmp(argv[i],"--vcf") == 0 || strcmp(argv[i],"-v") == 0 ) {
            focal_population_vcf_file = argv[++i] ;
        }
        if ( strcmp(argv[i],"--vcf") == 0 || strcmp(argv[i],"-v") == 0 ) {
            focal_population_vcf_file = argv[++i] ;
        }
        if ( strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--cds") == 0) {
            cds_file = argv[++i] ;
        }
        if ( strcmp(argv[i],"-i") == 0 ) {
            output_invariant = false ;
        }
    }
    if ( focal_population_vcf_file == "null" || outgroup_file == "null" ) {
        cerr << "MUST SPECIFY VCF AND AT LEAST ONE OUTGROUP FILE\n\n\n" ;
        exit(1) ;
    }
    if ( cds_file == "null" ) {
        cerr << "MUST SPECIFY ANNOTATION FILE\n\n\n" ;
        exit(1) ;
    }
}

//// functions for splitting info fields
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

//////// CDS FUNCTIONS AND INFORMATION

class poly_div_data {
public:
    int position ;
    string reference_allele ;
    string alternate_allele ;
    int alt_counts ;
    int number_chromosomes ;
    string outgroup_allele ;
    int outgroup_allele_score ;
    
} ;

class cds_entry {
    
public:
    
    //// CDS INFORMATION
    string cds_id ;
    string contig ;
    int start ;
    int stop ;
    int offset ;
    string strand ;
    
    /// will hold id's for parent transcripts to be converted into mk tables later via summation 
    vector<string> parent_transcripts ;
    
    //// FUNCTION TO EXTRACT CDS INFROMATION
    void parse_cds_line ( ifstream &CDS ) ;
    ///  FUNCTION TO SORT CDS INFORMATION
    friend bool operator<(cds_entry a, cds_entry b) {
        return a.start < b.start ;
    }
    
    //// polymorphism and divergence information
    int outgroup_iterator ;
    vector<poly_div_data> sequence_data ;
    
} ;

void cds_entry::parse_cds_line ( ifstream &CDS ) {
    
    string trash ;                  /// placeholder for information that we don't need
    
    CDS >> contig ;
    CDS >> trash ;
    CDS >> trash ;
    CDS >> start ;
    CDS >> stop ;
    CDS >> trash ;
    CDS >> strand ;
    CDS >> offset ;

    getline(CDS,cds_id) ;

}

vector<string> identify_parent_transcripts( string &cds_file, string &cds_info ) {

    vector<string> info_fields = split( cds_info, ';' ) ;
    vector<string> parent_ids ;

    //// gff3 or gff file
    if ( cds_file.substr(cds_file.size()-3) == "gff" || cds_file.substr(cds_file.size()-4) == "gff3" ) {
        /// should have fields: Parent
        for ( int f = 0 ; f < info_fields.size() ; f ++ ) {
            if ( info_fields.at(f).substr(0,7) == "Parent=" ) {
                parent_ids = split ( info_fields.at(f).substr(7), ',' ) ;
            }
            else if ( info_fields.at(f).substr(1,7) == "Parent=" ) {
                parent_ids = split ( info_fields.at(f).substr(8), ',' ) ;
            }
        }
    }

    //// gtf file
    else if ( cds_file.substr(cds_file.size()-3) == "gtf" ) {
        for ( int f = 0 ; f < info_fields.size() ; f ++ ) {
            if ( info_fields.at(f).substr(0,15) == " transcript_id " ) {
                parent_ids.push_back( info_fields.at(f).substr(16,info_fields.at(f).size()-17) ) ;
            }
        }
    }
    return parent_ids ;
}

void read_cds_file( string &cds_file, map< string, vector<cds_entry> > &cds_information ) {
    
    ifstream cds_stream ( cds_file.c_str() ) ;
    string cds_line ;
    while( !cds_stream.eof() ) {
        cds_entry new_cds ;
        new_cds.parse_cds_line( cds_stream ) ;
        new_cds.parent_transcripts = identify_parent_transcripts( cds_file, new_cds.cds_id ) ;
        new_cds.outgroup_iterator = 0 ;
        cds_information[new_cds.contig].push_back(new_cds) ;
    }
}

void collapse_redudant_cds ( map< string, vector<cds_entry> > &cds_information ) {
    
    /// for iterating through paired contig names / cds entries
    typedef map<string, vector<cds_entry> >::iterator it_type ;
    
    /// iterate through backwards and erase redundant cds copies
    for( it_type iterator = cds_information.begin(); iterator != cds_information.end(); iterator++ ) {
        for ( int i = iterator->second.size() - 1 ; i > 0 ; i -- ) {
           if ( iterator->second.at(i).start == iterator->second.at(i-1).start &&
                iterator->second.at(i).stop == iterator->second.at(i-1).stop &&
                iterator->second.at(i).strand == iterator->second.at(i-1).strand &&
                iterator->second.at(i).offset == iterator->second.at(i-1).offset ) {

                //// record additional parents
                for ( int p = 0 ; p < iterator->second.at(i).parent_transcripts.size() ; p ++ ) {
                    iterator->second.at(i-1).parent_transcripts.push_back( iterator->second.at(i).parent_transcripts.at(p) ) ;
                }
                
                /// now remove the second entry
                iterator->second.erase( iterator->second.begin() + i ) ;
            }
        }
    }
}

/// sort all cds entries by start coordinates vcf's go through in order so this is easy
void sort_cds_information( map< string, vector<cds_entry> > &cds_information ) {
    
    /// for iterating through paired contig names / cds entries
    typedef map<string, vector<cds_entry> >::iterator it_type ;
    
    //// iterate through each contig and sort cds by start position
    for( it_type iterator = cds_information.begin(); iterator != cds_information.end(); iterator++ ) {
        sort( iterator->second.begin(), iterator->second.end() ) ;
    }
}


/// read vcf data
class vcf_entry {
    
public:
    
    string contig ;
    int position ;                  /// one-based in vcf format
    string ref ;
    string alt ;
    float qual ;
    string info ;
    string format ;
    
    void read_vcf_line ( ifstream &VCF ) ;
    
} ;

void vcf_entry::read_vcf_line ( ifstream &VCF ) {

    string vcf_line ;                  /// placeholder for information that we don't need
    
    getline(VCF, vcf_line) ; 
    istringstream vcf_stream ( vcf_line ) ; 

    vcf_stream >> contig ;
    /// skip comment lines, and skip empty contigs
    if ( contig.size() == 0 || contig.at(0) == '#' ) {
        return ;
    }
    
    string trash ; 

    vcf_stream >> position ;
    vcf_stream >> trash ;
    vcf_stream >> ref ;
    vcf_stream >> alt ;
    vcf_stream >> trash ;
    qual = atoi( trash.c_str() ) ; 
    vcf_stream >> trash ;
    vcf_stream >> info ;

}

int extract_AN ( vector<string> &info ) {
    for ( int v = 0 ; v < info.size() ; v ++ ) {
        if ( info.at(v).substr(0,3) == "AN=" ) {
            return atoi( info.at(v).substr(3).c_str() ) ;
        }
    }
    //// if we don't find AN, return 0.
    //// maybe should throw a fail, since this should never happen in a VCF file
    return 0 ;
}

vector<int> extract_AC ( vector<string> &info ) {
    for ( int v = 0 ; v < info.size() ; v ++ ) {
        if ( info.at(v).substr(0,3) == "AC=" ) {
            vector<string> AF_string = split( info.at(v).substr(3), ',' ) ;
            vector<int> AF ;
            for ( int i = 0 ; i < AF_string.size() ; i ++ ) {
                AF.push_back(atoi(AF_string.at(i).c_str())) ;
            }
            return AF ;
        }
    }
    vector<int> fail ;
    fail.push_back(0) ;
    return fail ;
}

/// parse relevant vcf information for polymorphism dataset
poly_div_data parse_vcf_information ( vcf_entry &new_vcf ) {
    
    poly_div_data new_site ;
    new_site.position = new_vcf.position ;
    new_site.reference_allele = new_vcf.ref ;
    
    /// parse relevant vcf information
    vector<string> vcf_line_info = split( new_vcf.info, ';' ) ;

    /// extract AN
    new_site.number_chromosomes = extract_AN ( vcf_line_info ) ;
    
    //// if the site is invariant reference, just this info is needed
    if ( new_vcf.alt == "." ) {
        new_site.alt_counts = 0 ;
        return new_site ;
    }
    
    //// otherwise, we will parse relevant information about variable sites from VCF line
    vector<int> counts = extract_AC( vcf_line_info ) ;
    
    /// if there is only one alt allele freq, we will just record that
    if ( counts.size() == 1 && new_site.number_chromosomes > counts.at(0) ) {
        new_site.alt_counts = counts.at(0) ;
        new_site.alternate_allele = new_vcf.alt ;
        return new_site ;
    }

    //// unless there is only one alternate allele and no reference allele, then we simply swap that 
    if ( new_site.number_chromosomes == counts.at(0) ) { 
        new_site.alt_counts = 0 ; 
        new_site.reference_allele = new_vcf.alt ; 
        return new_site ; 
    }    
    
    //// extract alternate alleles
    vector<string> alternate_alleles = split( new_vcf.alt, ',' ) ;
    
    ///  the reference allele is not present
    if ( counts.at(0) + counts.at(1) == new_site.number_chromosomes ) {
        new_site.alternate_allele = alternate_alleles.at(1) ;
        new_site.alt_counts = counts.at(1) ;
        new_site.reference_allele = alternate_alleles.at(0) ;
        return new_site ;
    }
    
    /// if triallelic, we will retain the two more frequent sites
    int reference_count = new_site.number_chromosomes - counts.at(0) - counts.at(1) ;
    if ( counts.at(1) < reference_count && counts.at(1) < counts.at(0) ) {
        new_site.alternate_allele = alternate_alleles.at(0) ;
        new_site.alt_counts = counts.at(0) ;
        new_site.number_chromosomes -= counts.at(1) ;
        return new_site ;
    }
    if ( counts.at(0) < reference_count && counts.at(0) < counts.at(1) ) {
        new_site.alternate_allele = alternate_alleles.at(1) ;
        new_site.alt_counts = counts.at(1) ;
        new_site.number_chromosomes -= counts.at(0) ;
        return new_site ;
    }
    else {
        new_site.alternate_allele = alternate_alleles.at(1) ;
        new_site.reference_allele = alternate_alleles.at(0) ;
        new_site.number_chromosomes = counts.at(0) + counts.at(1) ;
        new_site.alt_counts = counts.at(1) ;
        return new_site ;
    }
}

void import_vcf_data ( string vcf_file, map< string, vector<cds_entry> > &cds_information, float &min_qual ) {
    
    int cds_index ;
    string contig = "null" ;
    
    ifstream VCF ( vcf_file.c_str() ) ;
    while( !VCF.eof() ) {
        
        /// import new vcf line
        vcf_entry new_vcf ;
        new_vcf.read_vcf_line( VCF ) ;

        /// skip comment lines
        if ( new_vcf.contig.size() == 0 || new_vcf.contig.at(0) == 0 ) {
            continue ;
        }
        
        /// If vcf entry is below the minimum quality, skip to the next one
        if ( new_vcf.qual < min_qual ) {
            continue ; 
        }
        
        /// if vcf has skipped to the next chromosome, reset the cds index to zero
        if ( contig != new_vcf.contig && new_vcf.contig.size() != 0 ) {
            cerr << "\t\t\tcontig=" << new_vcf.contig << endl ; 
            contig = new_vcf.contig ;
            cds_index = 0 ;
        }
        
        /// if the vcf is past the stop of the current cds, skip to next cds
        while ( cds_index < cds_information[contig].size() && new_vcf.position > cds_information[contig].at(cds_index).stop ) {
            cds_index ++ ;
        }
        
        /// skip if we've scrolled past the last cds until we get to the next contig
        if ( cds_index == cds_information[contig].size() ) {
            continue ;
        }
        
        /// if the vcf is behind the cds just go to the next vcf entry
        if ( new_vcf.position < cds_information[contig].at(cds_index).start ) {
            continue ;
        }
        
        /// if those conditions are met, the vcf line is somewhere within the range of the current cds (and possibly w/in the next few cds)
        /// first, parse the relevant vcf information
        poly_div_data new_site = parse_vcf_information( new_vcf ) ;
        new_site.outgroup_allele = "N" ;
        new_site.outgroup_allele_score = 0 ;    /// this should only matter for maf outgroup files
        
        //// set end of cds checking
        int stop = 10 ;             /// be default we'll check the next 10 cds
        if ( cds_index + stop >= cds_information[contig].size() ) {
            stop = cds_information[contig].size() - cds_index - 1 ;
        }
        
        /// check to ensure we are within the coordinates of this CDS ( will always be true when c == 0, but for later CDS's this needs to exist )
        for ( int c = 0 ; c < stop ; c ++ ) {
            if ( new_vcf.position >= cds_information[contig].at(cds_index+c).start && new_vcf.position <= cds_information[contig].at(cds_index+c).stop ) {
                cds_information[contig].at(cds_index+c).sequence_data.push_back(new_site) ;
            }
        }
    }
}

//// identify outgroup alleles
vector<string> outgroup_vcf_information( vcf_entry &new_vcf ) {
    
    vector<string> outgroup_alleles ;
    
    //// if the site is invariant reference, just this info is needed
    if ( new_vcf.alt == "." ) {
        outgroup_alleles.push_back( new_vcf.ref ) ;
        return outgroup_alleles ;
    }
    
    /// okay if there's an alterate allele
    vector<string> vcf_line_info = split( new_vcf.info, ';' ) ;
    int AN = extract_AN( vcf_line_info ) ;
    vector<int> AC = extract_AC( vcf_line_info ) ;
    vector<string> alternate_alleles = split( new_vcf.alt, ',' ) ;
    
    /// if there are two alternates
    if ( AC.size() == 2 ) {
        return alternate_alleles ;
    }
    /// if there is one fixed difference
    if ( AN == AC.at(0) ) {
        return alternate_alleles ;
    }
    /// otherwise, this is a heterozygous site
    alternate_alleles.push_back( new_vcf.ref ) ;
    return alternate_alleles ;
}

//// okay, now import information from outgroup VCF file
void import_outgroup_vcf ( string vcf_file, map< string, vector<cds_entry> > &cds_information, float &min_qual ) {
    
    int cds_index ;
    string contig = "null" ;
    
    ifstream VCF ( vcf_file.c_str() ) ;

    while( !VCF.eof() ) {
        
        /// read vcf line
        vcf_entry new_vcf ;
        new_vcf.read_vcf_line( VCF ) ;
        
        /// skip comment line
        if ( new_vcf.contig.size() == 0 || new_vcf.contig.at(0) == '#' ) {
            continue ;
        }
        
        ///// skip if below minimum quality
        if ( new_vcf.qual < min_qual ) {
            continue ;
        }
        
        /// if vcf has skipped to the next chromosome, reset the cds index to zero
        if ( contig != new_vcf.contig && new_vcf.contig.size() != 0 ) {
            cerr << "\t\t\tcontig=" << new_vcf.contig << endl ; 
            contig = new_vcf.contig ;
            cds_index = 0 ;
        }

        /// if the vcf is past the stop of the current cds, skip to next cds
        while ( cds_index < cds_information[contig].size() && new_vcf.position > cds_information[contig].at(cds_index).stop ) {
            cds_index ++ ;
        }

        /// skip if we've scrolled past the last cds until we get to the next contig
        if ( cds_index == cds_information[contig].size() ) {
            continue ;
        }
        
        /// if the vcf is behind the cds just go to the next vcf entry
        if ( new_vcf.position < cds_information[contig].at(cds_index).start ) {
            continue ;
        }
        
        /// okay, so again, we're in the money here, now we need to find what the outgroup base is
        vector<string> outgroup_bases = outgroup_vcf_information( new_vcf ) ;
        
        //// set end of cds checking
        int stop = 10 ;             /// be default we'll check the next 10 cds
        if ( cds_index + stop >= cds_information[contig].size() ) {
            stop = cds_information[contig].size() - cds_index - 1 ;
        }

        /// okay, now we need to place this in the cds file that overlap this position
        for ( int c = 0 ; c < stop ; c ++ ) {

            ///// skip CDSs with fewer than 3 bp in focal population
            if ( cds_information[contig].at(cds_index+c).sequence_data.size() < 3 ) {
                continue ;
            }

            if ( new_vcf.position >= cds_information[contig].at(cds_index+c).start && new_vcf.position <= cds_information[contig].at(cds_index+c).stop ) {
                
                //// iterate up until we're at the same position for ingroup and outgroup alleles
                while ( cds_information[contig].at(cds_index+c).outgroup_iterator < cds_information[contig].at(cds_index+c).sequence_data.size() - 1 && cds_information[contig].at(cds_index+c).sequence_data.at( cds_information[contig].at(cds_index+c).outgroup_iterator ).position < new_vcf.position ) {
                    cds_information[contig].at(cds_index+c).outgroup_iterator ++ ;
                }

                if ( cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).position == new_vcf.position ) {
                    
                    //// if there is only one outgroup allele, that's it
                    if ( outgroup_bases.size() == 1 ) {
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).outgroup_allele = outgroup_bases.at(0) ;
                        continue ;
                    }
                    
                    /// if there are two outgroup alleles and the first matches the reference base, assign it
                    if ( outgroup_bases.at(0) == cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele ) {
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).outgroup_allele = outgroup_bases.at(0) ;
                        continue ;
                    }
                    
                    //// switch reference base if it matches the outgroup base
                    /// i.e. reference is redefined to ancestral.
                    if ( outgroup_bases.at(0) == cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alternate_allele ) {
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).outgroup_allele = outgroup_bases.at(0) ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alternate_allele = cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele = outgroup_bases.at(0) ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alt_counts = cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).number_chromosomes - cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alt_counts ;
                        continue ;
                    }
                    
                    if ( outgroup_bases.at(1) == cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele ) {
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).outgroup_allele = outgroup_bases.at(1) ;
                        continue ;
                    }
                    
                    //// again switch the reference allele to match the outgroup
                    if ( outgroup_bases.at(1) == cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alternate_allele ) {
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).outgroup_allele = outgroup_bases.at(1) ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alternate_allele = cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele = outgroup_bases.at(1) ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alt_counts = cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).number_chromosomes - cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alt_counts ;
                        continue ;
                    }
                    
                    //// finally if neither allele is matched... we will just take the first one
                    cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).outgroup_allele = outgroup_bases.at(0) ;
                    
                    /// and make the reference (ancestral) site the more frequent of the two polymorphic sites
                    if ( cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alt_counts * 2 > cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).number_chromosomes ) {
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alt_counts = cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).number_chromosomes - cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alt_counts ;
                        string holder = cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alternate_allele ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).alternate_allele = cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele ;
                        cds_information[contig].at(cds_index+c).sequence_data.at(cds_information[contig].at(cds_index+c).outgroup_iterator).reference_allele = holder ;
                    }
                }
            }
        }
    }
}

class maf_entry {
  
public:
    
    string contig ;
    int start ;
    int alignment_length ;
    string reference_strand ;
    string reference_sequence ;
    string outgroup_strand ;
    string outgroup_sequence ;
    
    void read_maf_line( ifstream &MAF ) ;
    
} ;

/// reads next two lines ( i.e. we assume a single reference and single query, other maf's are not accepted )
void maf_entry::read_maf_line ( ifstream &MAF ) {

    string trash ;
    
    /// import reference line
    MAF >> trash >> contig >> start >> alignment_length >> reference_strand ;
    MAF >> trash >> reference_sequence ;

    /// import outgroup line
    MAF >> trash ;
    MAF >> trash ;
    MAF >> trash ;
    MAF >> trash >> outgroup_strand ;
    MAF >> trash >> outgroup_sequence ;
    
    /// convert to 1-based coordinates since gff and vcf are 1-based
    start ++ ;

}

/// SLIGHTLY DIFFERENT FUNCTION TO IMPORT MAF FORMATTED OUTGROUP FILE
void import_outgroup_maf ( string maf_file, map< string, vector<cds_entry> > &cds_information, float &min_qual, bool &accept_lower ) {
    
    ifstream MAF ( maf_file.c_str() ) ;
    string maf_line ;
    while( getline( MAF, maf_line ) ) {
        
        //// this is the first line in an alignment block
        if ( maf_line.size() > 0 && maf_line.at(0) == 'a' ) {
            
            //// get the alignment score
            float score = atof( maf_line.substr(8).c_str() ) ;
            
            //// check if the alignment score is below our minimum threshhold
            if ( score < min_qual ) {
                continue ;
            }
            
            maf_entry new_maf ;
            new_maf.read_maf_line( MAF ) ;

            //// check if this contig exists in the set of cdses
            if ( cds_information.count(new_maf.contig) == 0 ) {
                continue ;
            }

            //// check all cds entries on the reference contig for overlap
            /// THIS IS A SLOW IMPLEMENTATION, WE COULD TRY TO PRESORT THE MAF FILE
            for ( int c = 0 ; c < cds_information[new_maf.contig].size() ; c ++ ) {
                
                //// skip those cases with no focal population polymorphisms
                if ( cds_information[new_maf.contig].at(c).sequence_data.size() == 0 ) {
                    continue ;
                }

                //// check for overlap between ranges
                if ( ( new_maf.start >= cds_information[new_maf.contig].at(c).start && new_maf.start <= cds_information[new_maf.contig].at(c).stop ) || ( new_maf.start + new_maf.alignment_length >= cds_information[new_maf.contig].at(c).start && new_maf.start + new_maf.alignment_length <= cds_information[new_maf.contig].at(c).stop ) || ( new_maf.start <= cds_information[new_maf.contig].at(c).start && new_maf.start + new_maf.alignment_length >= cds_information[new_maf.contig].at(c).stop ) ) {
                    
                    /// now scroll through alignment and add outgroup sites to focal population sites
                    int alignment_index = -1 ;
                    int cds_index = 0 ;
                    for ( int s = 0 ; s < new_maf.outgroup_sequence.size() ; s ++ ) {
                   
                        //// exclude those where the ingroup is an indel
                        if ( new_maf.reference_sequence.at(s) == '-' ) {
                            continue ;
                        }
                        /// augment site count if not an indel in the ingroup
                        alignment_index ++ ;

                        /// skip sites where the outgroup is an indel
                        if ( new_maf.outgroup_sequence.at(s) == '-' ) {
                            continue ;
                        }
                   
                        /// skip sites wherein we havent' caught up to the position in the cds
                        if ( alignment_index + new_maf.start < cds_information[new_maf.contig].at(c).sequence_data.at(cds_index).position ) {
                            continue ;
                        }

                        /// if we're past the current cds, scroll the next one
                        while ( cds_index < cds_information[new_maf.contig].at(c).sequence_data.size() && alignment_index + new_maf.start > cds_information[new_maf.contig].at(c).sequence_data.at(cds_index).position ) {
                            cds_index ++ ;
                        }
                    
                        /// break out of the loop if the alignment went past the last position within the cds
                        if ( cds_index == cds_information[new_maf.contig].at(c).sequence_data.size() ) {
                            break ;
                        }
                    
                        /// if upper-case is required, skip soft-masked bases
                        if ( !isupper( new_maf.outgroup_sequence.at(s) ) && accept_lower == false ) { 
                            continue ; 
                        }

                        /// okay, now that we're at the same position
                            /// only record upper-case bases (i.e. non-softmasked in most assemblies)
                            /// only add if this base has the higher score at this site
                        if ( cds_information[new_maf.contig].at(c).sequence_data.at(cds_index).position == alignment_index + new_maf.start && cds_information[new_maf.contig].at(c).sequence_data.at(cds_index).outgroup_allele_score < score ) {
                            cds_information[new_maf.contig].at(c).sequence_data.at(cds_index).outgroup_allele  = toupper( new_maf.outgroup_sequence.at(s) ) ;
                            cds_information[new_maf.contig].at(c).sequence_data.at(cds_index).outgroup_allele_score = score ;
                        }
                    }
                }
            }
        }
    }
}

//// return true for maf
//// return false for vcf
bool check_file_type( string &file ) {
    if ( file.substr(file.size()-3) == "maf" ) {
        return 1 ;
    }
    return 0 ;
}

//// compliment sequence
string compliment_sequence ( const string &plus_strand_sequence ) {
    string minus_strand_sequence ;
    for ( int l = 0 ; l < plus_strand_sequence.size() ; l ++ ) {
        if ( plus_strand_sequence.at(l) == 'A' ) {
            minus_strand_sequence += "T" ;
            continue ;
        }
        if ( plus_strand_sequence.at(l) == 'T' ) {
            minus_strand_sequence += "A" ;
            continue ;
        }
        if ( plus_strand_sequence.at(l) == 'G' ) {
            minus_strand_sequence += "C" ;
            continue ;
        }
        if ( plus_strand_sequence.at(l) == 'C' ) {
            minus_strand_sequence += "G" ;
        }
    }
    return minus_strand_sequence ;
}

//// will scroll across each cds entry and output stats at relevant sites
void compute_mk_info( map< string, vector<cds_entry> > &cds_information, codon_transition_counts &transition_counts, bool output_invariant ) {

    /// iterate across all chromosomes
    for ( auto chromosome = cds_information.begin(); chromosome != cds_information.end(); chromosome++ ) {
        
        /// iterate across all cds's
        for ( int c = 0 ; c < chromosome->second.size() ; c ++ ) {

            //// we must have at least three bases to have any meaningful sequence data
            if ( chromosome->second.at(c).sequence_data.size() < 3 ) {
                continue ;
            }
            
            /// on the + strand, we don't have to do much
            if ( chromosome->second.at(c).strand == "+" ) {
                
                /// iterate through the vector of stored polymorphism data
                int cds_position_iterator = 0 ;
                for ( int p = chromosome->second.at(c).start + chromosome->second.at(c).offset ; p < chromosome->second.at(c).stop - 1 ; p += 3 ) {

                    /// if we haven't caught up up to this position
                    while ( cds_position_iterator < chromosome->second.at(c).sequence_data.size() - 1 &&  chromosome->second.at(c).sequence_data.at(cds_position_iterator).position < p ) {
                        cds_position_iterator ++ ;
                    }

                    //// alternatively, we might be passed this site, then skip to the next codon
                    if ( chromosome->second.at(c).sequence_data.at(cds_position_iterator).position > p ) {
                        continue ;
                    }

                    //// check to make sure all bases in this codon have a polymorphic base
                    if ( cds_position_iterator + 2 >= chromosome->second.at(c).sequence_data.size() || chromosome->second.at(c).sequence_data.at(cds_position_iterator).position != p || chromosome->second.at(c).sequence_data.at(cds_position_iterator+1).position != p + 1 || chromosome->second.at(c).sequence_data.at(cds_position_iterator+2).position != p + 2 ) {
                        continue ;
                    }

                    /// check to make sure each site has an outgroup 1 allele called as well
                    if ( chromosome->second.at(c).sequence_data.at(cds_position_iterator).outgroup_allele == "N" || chromosome->second.at(c).sequence_data.at(cds_position_iterator+1).outgroup_allele == "N" || chromosome->second.at(c).sequence_data.at(cds_position_iterator+2).outgroup_allele == "N" ) {
                        continue ;
                    }
  
                    /// identify polymorphic sites and record dn and ds information
                    string base_codon ;
                    string outgroup_codon ;
                    for ( int s = cds_position_iterator ; s < cds_position_iterator + 3 ; s ++ ) {
                        base_codon += chromosome->second.at(c).sequence_data.at(s).reference_allele ;
                        outgroup_codon += chromosome->second.at(c).sequence_data.at(s).outgroup_allele ;
                    }
                    
                    float average_counts = 0 ;
                    
                    for ( int s = 0 ; s < 3 ; s ++ ) {
                        
                        /// compute polymorphic counts
                        int ps = 0 ;
                        int pn = 0 ;
                        string alternate_codon = base_codon ;
                        alternate_codon.replace(s,1,chromosome->second.at(c).sequence_data.at( cds_position_iterator + s ).alternate_allele ) ;
                        string transition = base_codon + alternate_codon ;
                        ps = transition_counts.synonymous[transition][s] ;
                        pn = transition_counts.non_synonymous[transition][s] ;

                        /// divergent counts
                        transition = base_codon + outgroup_codon ;
                        int dn = 0 ;
                        int ds = 0 ;
                        if ( outgroup_codon[s] != base_codon[s] && ( alternate_codon.size() < 3 || alternate_codon[s] != outgroup_codon[s] ) ) {
                            dn = transition_counts.non_synonymous[transition][s] ;
                            ds = transition_counts.synonymous[transition][s] ;
                        }
                        
                        /// if alternate codon equals outgroup switch counts
                        if ( alternate_codon[s] == outgroup_codon[s] ) {
                            alternate_codon[s] = base_codon[s] ;
                            base_codon[s] = outgroup_codon[s] ;
                            chromosome->second.at(c).sequence_data.at( cds_position_iterator + s ).alt_counts = chromosome->second.at(c).sequence_data.at( cds_position_iterator + s ).number_chromosomes - chromosome->second.at(c).sequence_data.at( cds_position_iterator + s ).alt_counts ;
                        }
                        
                        /// for all variable sites, print the relevant polymorphism/divergent information
                        if ( output_invariant == true || ( ds > 0 || dn > 0 || ps > 0 || pn > 0 ) ) {
                            cout << chromosome->first << "\t" << p + s << "\t" ;
                            cout << outgroup_codon[s] << "\t" << base_codon[s] << "\t" ;
                            if ( pn > 0 || ps > 0 ) {
                                cout << alternate_codon[s] << "\t" ;
                            }
                            else {
                                cout << "NA" << "\t" ;
                            }
                            cout << dn << "\t" << ds << "\t" << pn << "\t" << ps << "\t" ;
                            cout << chromosome->second.at(c).sequence_data.at( cds_position_iterator + s ).number_chromosomes << "\t" ;
                            if ( pn > 0 || ps > 0 ) {
                                cout << (float)chromosome->second.at(c).sequence_data.at( cds_position_iterator + s ).alt_counts/(float)chromosome->second.at(c).sequence_data.at( cds_position_iterator + s ).number_chromosomes << "\t" ;
                            }
                            else {
                                cout << 0 << "\t" ;
                            }
                            
                            //// now print , delimited list of all transcripts
                            for ( int t = 0 ; t < chromosome->second.at(c).parent_transcripts.size() ; t ++ ) {
                                cout << chromosome->second.at(c).parent_transcripts.at(t) << "," ;
                            }
                            cout << "\t" << fold_counts[base_codon][s] << endl ;
                        }
                    }
                }
            }
            
            //// now deal with the minus strand cds
            else {
                
                //// here we will iterate from back to front
                int cds_position_iterator = chromosome->second.at(c).sequence_data.size() - 1 ;
                
                for ( int p = chromosome->second.at(c).stop - chromosome->second.at(c).offset ; p > chromosome->second.at(c).start + 1 ; p -= 3 ) {
                    
                    /// if we haven't caught up up to this position
                    while ( cds_position_iterator > 0 && chromosome->second.at(c).sequence_data.at(cds_position_iterator).position > p ) {
                        cds_position_iterator -- ;
                    }
                    
                    //// alternatively, we might be passed this site, then skip to the next codon
                    if ( cds_position_iterator < 2 || chromosome->second.at(c).sequence_data.at(cds_position_iterator).position < p ) {
                        continue ;
                    }
                    
                    //// check to make sure all bases in this codon have a polymorphic base
                    if ( chromosome->second.at(c).sequence_data.at(cds_position_iterator).position != p || chromosome->second.at(c).sequence_data.at(cds_position_iterator-1).position != p - 1 || chromosome->second.at(c).sequence_data.at(cds_position_iterator-2).position != p - 2 ) {
                        continue ;
                    }
                    
                    /// check to make sure each site has an outgroup 1 allele called as well
                    if ( chromosome->second.at(c).sequence_data.at(cds_position_iterator).outgroup_allele == "N" || chromosome->second.at(c).sequence_data.at(cds_position_iterator-1).outgroup_allele == "N" || chromosome->second.at(c).sequence_data.at(cds_position_iterator-2).outgroup_allele == "N" ) {
                        continue ;
                    }
                    
                    /// identify polymorphic sites and record dn and ds information
                    string base_codon ;
                    string outgroup_codon ;
                    for ( int s = cds_position_iterator ; s > cds_position_iterator - 3 ; s -- ) {                        
                        base_codon += chromosome->second.at(c).sequence_data.at(s).reference_allele ;
                        outgroup_codon += chromosome->second.at(c).sequence_data.at(s).outgroup_allele ;
                    }
                    
                    for ( int s = 0 ; s < 3 ; s ++ ) {
                        
                        /// compute polymorphic counts
                        int ps = 0 ;
                        int pn = 0 ;
                        string alternate_codon = base_codon ;
                        alternate_codon.replace(s,1,chromosome->second.at(c).sequence_data.at( cds_position_iterator - s ).alternate_allele) ;
                        string transition = compliment_sequence(base_codon + alternate_codon) ;
                        ps = transition_counts.synonymous[transition][s] ;
                        pn = transition_counts.non_synonymous[transition][s] ;
                        
                        /// if alternate codon equals outgroup switch counts
                        if ( alternate_codon[s] == outgroup_codon[s] ) {
                            alternate_codon[s] = base_codon[s] ;
                            base_codon[s] = outgroup_codon[s] ;
                            chromosome->second.at(c).sequence_data.at( cds_position_iterator - s ).alt_counts = chromosome->second.at(c).sequence_data.at( cds_position_iterator - s ).number_chromosomes - chromosome->second.at(c).sequence_data.at( cds_position_iterator - s ).alt_counts ;
                        }
                        
                        //// now compute divergent sites
                        transition = compliment_sequence(base_codon + outgroup_codon) ;
                        int dn = 0 ;
                        int ds = 0 ;
                        if ( outgroup_codon[s] != base_codon[s] && ( alternate_codon.size() < 3 || alternate_codon[s] != outgroup_codon[s] ) ) {
                            dn = transition_counts.non_synonymous[transition][s] ;
                            ds = transition_counts.synonymous[transition][s] ;
                        }
                        
                        /// for all variable sites, print the relevant polymorphism/divergent information
                        if ( output_invariant == true || ( ds > 0 || dn > 0 || ps > 0 || pn > 0 ) ) {
                            cout << chromosome->first << "\t" << p - s << "\t" ;
                            cout << compliment_sequence( outgroup_codon.substr(s, 1).c_str() ) << "\t" << compliment_sequence( base_codon.substr(s, 1).c_str() ) << "\t" ;
                            if ( pn > 0 || ps > 0 ) {
                                cout << compliment_sequence( alternate_codon.substr(s, 1).c_str() ) << "\t" ;
                            }
                            else {
                                cout << "NA" << "\t" ;
                            }
                        
                            cout << dn << "\t" << ds << "\t" << pn << "\t" << ps << "\t" ;
                            cout << chromosome->second.at(c).sequence_data.at( cds_position_iterator - s ).number_chromosomes << "\t" ;
                            if ( pn > 0 || ps > 0 ) {
                                cout << (float)chromosome->second.at(c).sequence_data.at( cds_position_iterator - s ).alt_counts/(float)chromosome->second.at(c).sequence_data.at( cds_position_iterator - s ).number_chromosomes << "\t" ;
                            }
                            else {
                                cout << "NA" << "\t" ;
                            }
                            
                            //// now print , delimited list of all transcripts
                            for ( int t = 0 ; t < chromosome->second.at(c).parent_transcripts.size() ; t ++ ) {
                                cout << chromosome->second.at(c).parent_transcripts.at(t) << "," ;
                            }
                            cout << "\t" << fold_counts[compliment_sequence(base_codon)][s] << endl ;

                        }
                    }
                }
            }
        }
    }
}

int main ( int argc, char **argv ) {
	
    /// for time benchmarking
    clock_t t = clock();
    clock_t last_time = clock() ;
    
    //// read command line options
    cerr << "reading command line options\n" ;
    cmd_line options ;
    options.read_cmd_line( argc, argv ) ;
    
    //// create counts of mutations required for each transition all is stored in codon_counts.h
    //// this is a very simple flavor of dynamic programming
    cerr << "populating codon transition lists\n" ;
    codon_transition_counts transition_counts ;
    transition_counts.populate_non_syn() ;
    transition_counts.populate_syn() ;

    /// benchmarking
    cerr << "initialization runtime:\t" << (double) (clock() - last_time) / 1000000 << endl ;
    last_time = clock() ;
    
    //// read through cds file and store relevant information
    cerr << "acquiring annotation information\n\t\tFILE=" << options.cds_file << endl ;
    map< string, vector<cds_entry> > cds_information ;
    read_cds_file( options.cds_file, cds_information ) ;
    sort_cds_information( cds_information ) ;
    
    //// this will remove cds entries that are completely redudant, e.g. those used in several transcripts
    cerr << "collapsing redundant CDS entries\n" ;
    collapse_redudant_cds( cds_information ) ;
    
    /// benchmarking
    cerr << "annotation runtime:\t" << (double) (clock() - last_time) / 1000000 << endl ;
    last_time = clock() ;
    
    //// now import vcf data for focal population
    cerr << "importing vcf data from focal population\n\t\tFILE=" << options.focal_population_vcf_file << endl ; 
    import_vcf_data( options.focal_population_vcf_file, cds_information, options.minimum_focal_vcf_quality ) ;
    
    /// benchmarking
    cerr << "vcf import runtime:\t" << (double) (clock() - last_time) / 1000000 << endl ;
    last_time = clock() ;
    
    //// now import outgroup alignment data for the outgroup
    cerr << "importing outgroup data\n\t\tFILE=" << options.outgroup_file << endl ; 
    bool type = check_file_type( options.outgroup_file ) ;
    if ( type == 0 ) {
        import_outgroup_vcf( options.outgroup_file, cds_information, options.minimum_outgroup_quality ) ;
    }
    else {
        import_outgroup_maf( options.outgroup_file, cds_information, options.minimum_outgroup_quality, options.accept_lower ) ;
    }
    
    /// benchmarking
    cerr << "outgroup import runtime:\t" << (double) (clock() - last_time) / 1000000 << endl ;
    last_time = clock() ;
    
    /// compute and output polymorphism and divergence information
    cerr << "computing MK information\n" ;
    compute_mk_info ( cds_information, transition_counts, options.output_invariant ) ;
    
    /// benchmarking
    cerr << "mk runtime:\t" << (double) (clock() - last_time) / 1000000 << endl ;
    
    /// output time
    cerr << "total runtime:\t" << (double) (clock() - t) / 1000000 << endl ;
    
    return(0) ;
}
