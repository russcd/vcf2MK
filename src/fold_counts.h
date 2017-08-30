//
//  fold_counts.h
//

//// here we have stored the number of steps between any pair of codons
#include <utility>
#include <map>
#include <string>
#include <vector> 

using namespace std;

#ifndef _fold_counts_h
#define _fold_counts_h

map<string,map<int,int> > create_counts () {
    
    map<string,map<int,int> > fold_counts ;
    
	fold_counts["TTG"][0] = 2 ;
	fold_counts["TTG"][1] = 0 ;
	fold_counts["TTG"][2] = 2 ;
	fold_counts["GTA"][0] = 0 ;
	fold_counts["GTA"][1] = 0 ;
	fold_counts["GTA"][2] = 4 ;
	fold_counts["CCG"][0] = 0 ;
	fold_counts["CCG"][1] = 0 ;
	fold_counts["CCG"][2] = 4 ;
	fold_counts["CCA"][0] = 0 ;
	fold_counts["CCA"][1] = 0 ;
	fold_counts["CCA"][2] = 4 ;
	fold_counts["TAC"][0] = 0 ;
	fold_counts["TAC"][1] = 0 ;
	fold_counts["TAC"][2] = 2 ;
	fold_counts["CCC"][0] = 0 ;
	fold_counts["CCC"][1] = 0 ;
	fold_counts["CCC"][2] = 4 ;
	fold_counts["GTG"][0] = 0 ;
	fold_counts["GTG"][1] = 0 ;
	fold_counts["GTG"][2] = 4 ;
	fold_counts["ACA"][0] = 0 ;
	fold_counts["ACA"][1] = 0 ;
	fold_counts["ACA"][2] = 4 ;
	fold_counts["CAG"][0] = 0 ;
	fold_counts["CAG"][1] = 0 ;
	fold_counts["CAG"][2] = 2 ;
	fold_counts["ATG"][0] = 0 ;
	fold_counts["ATG"][1] = 0 ;
	fold_counts["ATG"][2] = 0 ;
	fold_counts["GAT"][0] = 0 ;
	fold_counts["GAT"][1] = 0 ;
	fold_counts["GAT"][2] = 2 ;
	fold_counts["GCA"][0] = 0 ;
	fold_counts["GCA"][1] = 0 ;
	fold_counts["GCA"][2] = 4 ;
	fold_counts["TCG"][0] = 0 ;
	fold_counts["TCG"][1] = 0 ;
	fold_counts["TCG"][2] = 4 ;
	fold_counts["ATA"][0] = 0 ;
	fold_counts["ATA"][1] = 0 ;
	fold_counts["ATA"][2] = 3 ;
	fold_counts["TCT"][0] = 0 ;
	fold_counts["TCT"][1] = 0 ;
	fold_counts["TCT"][2] = 4 ;
	fold_counts["CGT"][0] = 0 ;
	fold_counts["CGT"][1] = 0 ;
	fold_counts["CGT"][2] = 4 ;
	fold_counts["ACT"][0] = 0 ;
	fold_counts["ACT"][1] = 0 ;
	fold_counts["ACT"][2] = 4 ;
	fold_counts["CGG"][0] = 2 ;
	fold_counts["CGG"][1] = 0 ;
	fold_counts["CGG"][2] = 4 ;
	fold_counts["GGA"][0] = 0 ;
	fold_counts["GGA"][1] = 0 ;
	fold_counts["GGA"][2] = 4 ;
	fold_counts["CAA"][0] = 0 ;
	fold_counts["CAA"][1] = 0 ;
	fold_counts["CAA"][2] = 2 ;
	fold_counts["TTC"][0] = 0 ;
	fold_counts["TTC"][1] = 0 ;
	fold_counts["TTC"][2] = 2 ;
	fold_counts["AAC"][0] = 0 ;
	fold_counts["AAC"][1] = 0 ;
	fold_counts["AAC"][2] = 2 ;
	fold_counts["AGA"][0] = 2 ;
	fold_counts["AGA"][1] = 0 ;
	fold_counts["AGA"][2] = 2 ;
	fold_counts["CTA"][0] = 2 ;
	fold_counts["CTA"][1] = 0 ;
	fold_counts["CTA"][2] = 4 ;
	fold_counts["AAG"][0] = 0 ;
	fold_counts["AAG"][1] = 0 ;
	fold_counts["AAG"][2] = 2 ;
	fold_counts["CTC"][0] = 0 ;
	fold_counts["CTC"][1] = 0 ;
	fold_counts["CTC"][2] = 4 ;
	fold_counts["GGC"][0] = 0 ;
	fold_counts["GGC"][1] = 0 ;
	fold_counts["GGC"][2] = 4 ;
	fold_counts["TGG"][0] = 0 ;
	fold_counts["TGG"][1] = 0 ;
	fold_counts["TGG"][2] = 0 ;
	fold_counts["GCC"][0] = 0 ;
	fold_counts["GCC"][1] = 0 ;
	fold_counts["GCC"][2] = 4 ;
	fold_counts["TGC"][0] = 0 ;
	fold_counts["TGC"][1] = 0 ;
	fold_counts["TGC"][2] = 2 ;
	fold_counts["GCG"][0] = 0 ;
	fold_counts["GCG"][1] = 0 ;
	fold_counts["GCG"][2] = 4 ;
	fold_counts["TCC"][0] = 0 ;
	fold_counts["TCC"][1] = 0 ;
	fold_counts["TCC"][2] = 4 ;
	fold_counts["ATT"][0] = 0 ;
	fold_counts["ATT"][1] = 0 ;
	fold_counts["ATT"][2] = 3 ;
	fold_counts["CCT"][0] = 0 ;
	fold_counts["CCT"][1] = 0 ;
	fold_counts["CCT"][2] = 4 ;
	fold_counts["AAA"][0] = 0 ;
	fold_counts["AAA"][1] = 0 ;
	fold_counts["AAA"][2] = 2 ;
	fold_counts["ACG"][0] = 0 ;
	fold_counts["ACG"][1] = 0 ;
	fold_counts["ACG"][2] = 4 ;
	fold_counts["TAT"][0] = 0 ;
	fold_counts["TAT"][1] = 0 ;
	fold_counts["TAT"][2] = 2 ;
	fold_counts["ACC"][0] = 0 ;
	fold_counts["ACC"][1] = 0 ;
	fold_counts["ACC"][2] = 4 ;
	fold_counts["TTA"][0] = 2 ;
	fold_counts["TTA"][1] = 0 ;
	fold_counts["TTA"][2] = 2 ;
	fold_counts["CAT"][0] = 0 ;
	fold_counts["CAT"][1] = 0 ;
	fold_counts["CAT"][2] = 2 ;
	fold_counts["CTT"][0] = 0 ;
	fold_counts["CTT"][1] = 0 ;
	fold_counts["CTT"][2] = 4 ;
	fold_counts["GAA"][0] = 0 ;
	fold_counts["GAA"][1] = 0 ;
	fold_counts["GAA"][2] = 2 ;
	fold_counts["CAC"][0] = 0 ;
	fold_counts["CAC"][1] = 0 ;
	fold_counts["CAC"][2] = 2 ;
	fold_counts["GCT"][0] = 0 ;
	fold_counts["GCT"][1] = 0 ;
	fold_counts["GCT"][2] = 4 ;
	fold_counts["GGT"][0] = 0 ;
	fold_counts["GGT"][1] = 0 ;
	fold_counts["GGT"][2] = 4 ;
	fold_counts["AGT"][0] = 0 ;
	fold_counts["AGT"][1] = 0 ;
	fold_counts["AGT"][2] = 2 ;
	fold_counts["AGG"][0] = 2 ;
	fold_counts["AGG"][1] = 0 ;
	fold_counts["AGG"][2] = 2 ;
	fold_counts["TCA"][0] = 0 ;
	fold_counts["TCA"][1] = 0 ;
	fold_counts["TCA"][2] = 4 ;
	fold_counts["ATC"][0] = 0 ;
	fold_counts["ATC"][1] = 0 ;
	fold_counts["ATC"][2] = 3 ;
	fold_counts["CTG"][0] = 2 ;
	fold_counts["CTG"][1] = 0 ;
	fold_counts["CTG"][2] = 4 ;
	fold_counts["GAG"][0] = 0 ;
	fold_counts["GAG"][1] = 0 ;
	fold_counts["GAG"][2] = 2 ;
	fold_counts["TGT"][0] = 0 ;
	fold_counts["TGT"][1] = 0 ;
	fold_counts["TGT"][2] = 2 ;
	fold_counts["CGC"][0] = 0 ;
	fold_counts["CGC"][1] = 0 ;
	fold_counts["CGC"][2] = 4 ;
	fold_counts["CGA"][0] = 2 ;
	fold_counts["CGA"][1] = 0 ;
	fold_counts["CGA"][2] = 4 ;
	fold_counts["GTT"][0] = 0 ;
	fold_counts["GTT"][1] = 0 ;
	fold_counts["GTT"][2] = 4 ;
	fold_counts["GTC"][0] = 0 ;
	fold_counts["GTC"][1] = 0 ;
	fold_counts["GTC"][2] = 4 ;
	fold_counts["AGC"][0] = 0 ;
	fold_counts["AGC"][1] = 0 ;
	fold_counts["AGC"][2] = 2 ;
	fold_counts["TTT"][0] = 0 ;
	fold_counts["TTT"][1] = 0 ;
	fold_counts["TTT"][2] = 2 ;
	fold_counts["AAT"][0] = 0 ;
	fold_counts["AAT"][1] = 0 ;
	fold_counts["AAT"][2] = 2 ;
	fold_counts["GGG"][0] = 0 ;
	fold_counts["GGG"][1] = 0 ;
	fold_counts["GGG"][2] = 4 ;
	fold_counts["GAC"][0] = 0 ;
	fold_counts["GAC"][1] = 0 ;
	fold_counts["GAC"][2] = 2 ;

    return fold_counts ;
}

map<string,map<int,int> > fold_counts = create_counts() ;

#endif
