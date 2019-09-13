#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

using namespace std;

int a, c, t, g;
int aa, ac, at, ag;
int ca, cc, ct, cg;
int ta, tc, tt, tg;
int ga, gc, gt, gg;
char cur;
char last;
double lines;
double letters;
void aSeq(char);
void cSeq(char);
void tSeq(char);
void gSeq(char);
void initializer();
void output();
void math(string);
string gaus();
double avglen, variance, deviation;
double aProb, cProb, gProb, tProb;
double aaProb, acProb, atProb, agProb;
double caProb, ccProb, ctProb, cgProb;
double taProb, tcProb, ttProb, tgProb;
double gaProb, gcProb, gtProb, ggProb;

int main(int argc, char *argv[]){
	std::ifstream file;
	string fileName;
	std::cout << "What file is the DNA stored in? " << flush;
	while(true){
		file.close();
    	file.clear();
		getline(cin, fileName); //get the file name from terminal
		file.open(fileName.c_str()); //open the file
		if(file){
			break;
		}
		std::cout << "Invalid file name, please enter a valid name: " << flush;
	}

	initializer();

	if (file.is_open()) {
	    std::string line;
		last = 'x'; //holds the last character (X if its the first character of the line)
		lines = 0;
		variance = 0;
	    while (getline(file, line)) {
			++lines;
			for(int i = 0; i<line.length(); ++i){
				++letters;
				cur = char(tolower(line.at(i))); //make sure the char is in lower case
				//count how many of each letter and letter combination there is
				if(cur == 'a'){
					++a;
					aSeq(last);
				}else if(cur == 'c'){
					++c;
					cSeq(last);
				}else if(cur == 't'){
					++t;
					tSeq(last);
				}else if(cur == 'g'){
					++g;
					gSeq(last);
				}
				last = cur;
			}
			last = 'x';
	    }
	    file.close();
		math(fileName);
		output();
	}
}

void math(string fileName){
	avglen = letters/lines;
	aProb = a/letters;
	cProb = c/letters;
	gProb = g/letters;
	tProb = t/letters;
	double totSeq = aa+ac+ag+at+ca+cc+cg+ct+ga+gc+gg+gt+ta+tc+tg+tt;
	aaProb = aa/totSeq;
	acProb = ac/totSeq;
	agProb = ag/totSeq;
	atProb = at/totSeq;
	caProb = ca/totSeq;
	ccProb = cc/totSeq;
	cgProb = cg/totSeq;
	ctProb = ct/totSeq;
	gaProb = ga/totSeq;
	gcProb = gc/totSeq;
	ggProb = gg/totSeq;
	gtProb = gt/totSeq;
	taProb = ta/totSeq;
	tcProb = tc/totSeq;
	tgProb = tg/totSeq;
	ttProb = tt/totSeq;

	//reopen the file to count how long each string is to calculate variance since the first time we looked through the file we didnt have the mean calculated yet
	std::ifstream file(fileName.c_str());
	if(file.is_open()){
		std::string line;
		while (getline(file, line)) {
			variance += pow(line.length() - avglen, 2);
		}
		file.close();
	}
	variance = variance/lines;
	deviation = sqrt(variance);
}

void output(){
	std::ofstream tofile("maxmiller.out");
	tofile << "Max Miller\n2317020\nMaxMiller@chapman.edu\n" << std::endl;

	tofile << "Sum:\nA : " << a << std::endl;
	tofile << "C : " << c << std::endl;
	tofile << "G : " << g << std::endl;
	tofile << "T : " << t << std::endl;
	tofile << "Total Lines : " << lines << std::endl;
	tofile << "Total Characters : " << letters << std::endl;
	tofile << "\nMean : " << avglen << std::endl;
	tofile << "Variance : " << variance << std::endl;
	tofile << "Standard Deviation : " << deviation << std::endl;
	tofile << "\nProbabilities:\nA : " << aProb << std::endl;
	tofile << "C : " << cProb << std::endl;
	tofile << "G : " << gProb << std::endl;
	tofile << "T : " << tProb << std::endl;
	tofile << "AA : " << aaProb << std::endl;
	tofile << "AC : " << acProb << std::endl;
	tofile << "AG : " << agProb << std::endl;
	tofile << "AT : " << atProb << std::endl;
	tofile << "CA : " << caProb << std::endl;
	tofile << "CC : " << ccProb << std::endl;
	tofile << "CG : " << cgProb << std::endl;
	tofile << "CT : " << ctProb << std::endl;
	tofile << "GA : " << gaProb << std::endl;
	tofile << "GC : " << gcProb << std::endl;
	tofile << "GG : " << ggProb << std::endl;
	tofile << "GT : " << gtProb << std::endl;
	tofile << "TA : " << taProb << std::endl;
	tofile << "TC : " << tcProb << std::endl;
	tofile << "TG : " << tgProb << std::endl;
	tofile << "TT : " << ttProb << "\n" << std::endl;

	for(int j = 0; j < 1000; ++j){
		tofile << gaus() << std::endl;
	}

	tofile.close();
}

string gaus(){
	//creates random sequences of DNA with the same Probabilitiesas the original file
	double x = ((double) rand() / (RAND_MAX));
	double y = ((double) rand() / (RAND_MAX));
	double bm = sqrt(-2*log(x))*cos(2*M_PI*y);
	int len = int(deviation*bm+avglen);

	int aChance = int(aProb*100);
	int cChance = int(cProb*100) + aChance;
	int gChance = int(gProb*100) + cChance;
	string str = "";
	for(int l = 0; l < len; ++l){
		int chance = rand() % 100 + 1;
		if(chance <= aChance){
			str += "A";
		}else if(chance <= cChance){
			str += "C";
		}else if(chance <= gChance){
			str += "G";
		}else{
			str += "T";
		}
	}
	return str;
}

void aSeq(char last){
	if(last == 'a'){
		++aa;
	}else if(last == 'c'){
		++ca;
	}else if(last == 't'){
		++ta;
	}else if(last == 'g'){
		++ga;
	}
}
void cSeq(char last){
	if(last == 'a'){
		++ac;
	}else if(last == 'c'){
		++cc;
	}else if(last == 't'){
		++tc;
	}else if(last == 'g'){
		++gc;
	}
}
void tSeq(char last){
	if(last == 'a'){
		++at;
	}else if(last == 'c'){
		++ct;
	}else if(last == 't'){
		++tt;
	}else if(last == 'g'){
		++gt;
	}
}
void gSeq(char last){
	if(last == 'a'){
		++ag;
	}else if(last == 'c'){
		++cg;
	}else if(last == 't'){
		++tg;
	}else if(last == 'g'){
		++gg;
	}
}

void initializer(){
	a=0;
	c=0;
	t=0;
	g=0;
	aa=0;
	ac=0;
	at=0;
	ag=0;
	ca=0;
	cc=0;
	ct=0;
	cg=0;
	ta=0;
	tc=0;
	tt=0;
	tg=0;
	ga=0;
	gc=0;
	gt=0;
	gg=0;
}
