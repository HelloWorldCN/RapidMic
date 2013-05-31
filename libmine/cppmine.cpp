#include <IOSTREAM>
#include <cstdlib>
#include "mine.h"
#include "cppmine.h"
#include "../src/stringenc.h"
#include "../src/arg_parser.h"
#include <fstream>
#include <time.h>
#if defined(WIN32) || defined(_WIN32)
#   include <windows.h>
#else
#   include <sys/time.h>
#endif
#if defined(WIN32) || defined(_WIN32)
int
	gettimeofday(struct timeval *tp, void *tzp)
{
	time_t clock;
	struct tm tm;
	SYSTEMTIME wtm;
	GetLocalTime(&wtm);
	tm.tm_year     = wtm.wYear - 1900;
	tm.tm_mon     = wtm.wMonth - 1;
	tm.tm_mday     = wtm.wDay;
	tm.tm_hour     = wtm.wHour;
	tm.tm_min     = wtm.wMinute;
	tm.tm_sec     = wtm.wSecond;
	tm. tm_isdst    = -1;
	clock = mktime(&tm);
	tp->tv_sec = clock;
	tp->tv_usec = wtm.wMilliseconds * 1000;
	return (0);
}
#endif



namespace {
	const char * const Program_name = "RapidMic";
	const char * const program_name = "rapidMic";
	const char * const program_year = "2013";
	const char * invocation_name = 0;
	void show_help()
	{
        std::printf( "---RapidMic---\n"
			"\nUsage: %s [options]\n", invocation_name );
		std::printf( "\nOptions:\n"
			"-h, --help     display this help and exit;\n"
			"-V, --version     output version information and exit;\n"
			"-i <file>, --input=<file>     input filename;\n"
            "-l <label style>, --inputlabel=<label style>     if input csv file's first line is comlumn name and the first comlumn is row name, then <label style>=0 ,else if just only comlumn name <label style>=1, else if just only row name <label style>=2 ;\n"
			"-a <alpha value>, --alpha=<alpha>     the exponent in B(n) = n^alpha (default: 0.6.) alpha must be in (0, 1.0]\n"
			"-c <clumps value>, --clumps=<c>     determines how many more clumps there will be than columns in every partition. Default value is 15,c must be > 0;\n"
			"-o <file>, --output=<file>     output filename (default: mine_out.csv);\n"
            "-L <label style>, --outputlabel=<label style>     output csv file adopt number index of row as vaiable label ,then <label style>=0，else adopt input file's row name, then <label style>=1;\n"
			"-A <allPairs>, --allPairs     will cause MINE to compare all pairs of variables against each other; \n"
			"-b <var index>, --pairsBetween=<var index>     will compare each of the first i variables to each of the rest of the variables. Variables are indexed from 0;input variable <var index> must be in (0, number of variables in file)\n"
			"-m <var index>, --master=<var index>     variable <var index> vs. all <var index> must be in [0, number of variables in file);\n"
			"-p <var1 index> <var2 index>,--onePair=<var1 index> <var2 index>     variable <var1 index> vs variable <var2 index> <var1 index> and <var2 index> must be in [0, number of variables in file);\n"		
			);

		std::printf( "\nReport bugs to tdm_2010#swjtu.edu.cn(replace # with @)\n"
			"Home page: http://www.swjtu.edu.cn\n" );
	}


	void show_version()
	{
		std::printf( "%s %s\n", Program_name, "SWJTU" );
	}


	void show_error( const char * const msg, const int errcode = 0,
		const bool help = false )
	{
		if( msg && msg[0] )
		{
			std::fprintf( stderr, "%s: %s", program_name, msg );
			if( errcode > 0 )
				std::fprintf( stderr, ": %s", std::strerror( errcode ) );
			std::fprintf( stderr, "\n" );
		}
		if( help )
			std::fprintf( stderr, "Try '%s --help' for more information.\n",
			invocation_name );
	}


	void internal_error( const char * const msg )
	{
		std::fprintf( stderr, "%s: internal error: %s.\n", program_name, msg );
		std::exit( 3 );
	}


	const char * optname( const int code, const Arg_parser::Option options[] )
	{
		static char buf[2] = "?";

		if( code != 0 )
			for( int i = 0; options[i].code; ++i )
				if( code == options[i].code )
				{ if( options[i].name ) return options[i].name; else break; }
				if( code > 0 && code < 256 ) buf[0] = code; else buf[0] = '?';
				return buf;
	}

} // end namespace
using namespace std;




/*
MINE constructor.

alpha is the exponent in B(n) = n^alpha and must be in (0,1], and c
determines how many more clumps there will be than columns in every
partition. c = 15 meaning that when trying to draw Gx grid lines on
the x-axis, the algorithm will start with at most 15*Gx clumps. c must
be > 0.

The constructor throws an exception when the parameters are invalid.
*/
MINE::MINE(double alpha, double c)
{
	char *ret;
	param.alpha = alpha;
	param.c = c;
	m_Results=NULL;
	m_ResultsArrayLen=0;
	m_inputDataMatrix=NULL;
	m_inputDataRowNum=0;
	m_inputDataColNum=0;
	m_masterId=m_betweenId=m_onepair1=m_onepair2=-1;
    m_outputLabel=m_inputLabel=-1;
	if ((ret=check_parameter(&param)))
		throw ret;
}


/*
MINE destructor.
*/
MINE::~MINE()
{
	if (m_Results!=NULL)
	{
		delete []m_Results;
	}

}


/*
Computes the maximum normalized mutual information scores between
the variables x and y of length n.
*/



/* Returns the Maximal Information Coefficient (MIC). */
double MINE::get_mic()
{
	if (m_Results==NULL&&m_AnalysisStyle!=OneParis)
	{
		return -1;
	}
	return m_Results->mic;
}


/* Returns the Maximum Asymmetry Score (MAS). */
double MINE::get_mas()
{
	if (m_Results==NULL&&m_AnalysisStyle!=OneParis)
	{
		return -1;
	}
	return m_Results->mas;
}


/* Returns the Maximum Edge Value (MEV). */
double MINE::get_mev()
{
	if (m_Results==NULL&&m_AnalysisStyle!=OneParis)
	{
		return -1;
	}
	return m_Results->mev;
}


/* Returns the Minimum Cell Number (MCN). */
double MINE::get_mcn()
{
	if (m_Results==NULL&&m_AnalysisStyle!=OneParis)
	{
		return -1;
	}
	return m_Results->mcn;
}

int MINE::AllPairsAnalysis( double **inData,int m,int n )
{
	m_AnalysisStyle=AllParis;
	if (m_Results!=NULL)
	{
		delete []m_Results;
	}
	int outlen=(m*(m-1))/2;
	m_Results=new mine_result_score[outlen];
	int ret=mine_allPairs_analysis(&param,inData,m,n,m_Results,outlen);
	if (ret)
	{
		m_inputDataRowNum=m;
		m_inputDataColNum=n;
		m_ResultsArrayLen=outlen;
	}else{//clear
		delete []m_Results;
		m_Results=NULL;
		m_ResultsArrayLen=0;

	}
	return ret;
	//cout <<m_Results[outlen-1].mic;

}
int MINE::TwoSetsAnalysis( double **inDataSet,int m,int n,int betweenid )
{
	m_AnalysisStyle=TwoSets;
	if (m_Results!=NULL)
	{
		delete []m_Results;
	}
	int outlen=betweenid*(m-betweenid);
	m_Results=new mine_result_score[outlen];
	int ret=mine_twoSetsAnalysis(&param,inDataSet, m, n,betweenid,m_Results,outlen);
	if (ret)
	{
		m_inputDataRowNum=m;
		m_inputDataColNum=n;
		m_ResultsArrayLen=outlen;
	}else{//clear
		delete []m_Results;
		m_Results=NULL;
		m_ResultsArrayLen=0;

	}

	return ret;
}
int MINE::MasterAnalysis( double **inData,int m,int n,int masterid )
{
	m_AnalysisStyle=MasterVariable;
	if (m_Results!=NULL)
	{
		delete []m_Results;
	}
	int outlen=m;
	m_Results=new mine_result_score[outlen];
	int ret=mine_masterVariableAnalysis(&param, inData, m, n, masterid, m_Results, outlen);
	if (ret)
	{
		m_inputDataRowNum=m;
		m_inputDataColNum=n;
		m_ResultsArrayLen=outlen;
	}else{//clear
		delete []m_Results;
		m_Results=NULL;
		m_ResultsArrayLen=0;

	}
	return ret;
}

int MINE::OnePairsAnalysis( double *x, double *y,int n )
{
	m_AnalysisStyle=OneParis;
	if (m_Results!=NULL)
	{
		delete []m_Results;
	}
	int outlen=1;
	m_Results=new mine_result_score[outlen];
	int ret=mine_onePairs_analysis(&param,x,y,n,m_Results);
	if (ret)
	{
		m_ResultsArrayLen=outlen;
	}else{//clear
		delete []m_Results;
		m_Results=NULL;
		m_ResultsArrayLen=0;

	}
	return ret;

}

int MINE::run( int argc, char **argv )
{
	if (!parserArgs(argc,argv)) return 0;
    //read input data file;
    cout << "begin reading input file;"<<endl;
    switch (m_inputLabel) {
        case 0:
            if(!readCSV(true, true)) return 0;
            break;
        case 1:
            if(!readCSV(false, true)) return 0;
            break;
        case 2:
            if(!readCSV(true, false)) return 0;
            break;
        default:
            break;
    }
    cout << "begin calculating......."<<endl;
    struct  timeval start;
    struct  timeval end;
    
    unsigned  long diff;
    gettimeofday(&start,NULL);
   
   
    int ret;
	switch (m_AnalysisStyle) {
	case AllParis:
		ret= AllPairsAnalysis(m_inputDataMatrix,m_inputDataRowNum,m_inputDataColNum);
		break;
	case TwoSets:
		if (m_betweenId>0&&m_betweenId<m_inputDataRowNum)
		{
			ret= TwoSetsAnalysis(m_inputDataMatrix,m_inputDataRowNum,m_inputDataColNum,m_betweenId);
		}else{ cout<<"input variable <var index> must be in (0, number of variables in file)"<<endl;return 0;}
		break;
	case MasterVariable:
		if (m_masterId>=0&&m_masterId<m_inputDataRowNum)
		{
			ret= MasterAnalysis(m_inputDataMatrix,m_inputDataRowNum,m_inputDataColNum,m_masterId);
		}else { cout<<"input variable <var index> must be in [0, number of variables in file)"<<endl;return 0;}
		break;
	case OneParis:
		if (m_onepair1<m_inputDataRowNum&&m_onepair1>=0&&m_onepair2>=0&&m_onepair2<m_inputDataRowNum)
		{
			ret= OnePairsAnalysis(m_inputDataMatrix[m_onepair1],m_inputDataMatrix[m_onepair2],m_inputDataColNum);
		}else { cout<<"input variable <var index> must be in [0, number of variables in file)"<<endl;return 0;}
	default:
		return 0;
		break;
	}
    if (ret) {
        cout << "completed calculation;"<<endl;
        gettimeofday(&end,NULL);
        diff = 1000000 * (end.tv_sec-start.tv_sec)+ end.tv_usec-start.tv_usec;
        cout <<"total calculating time:"<<diff/(1000*1000) <<"(s)"<<endl;
        return exportResult();
    }else return 0;
}

bool MINE::readCSV(bool rowlabel,bool collabel)
{
	vector<vector<double>> inputdata;
	releaseInputMatrix();
    m_inputDataRowNum=m_inputDataColNum=0;
    m_varColNames.clear();
    m_varRowNames.clear();
	if (!m_inputFileName.empty())
	{
		string line;
		ifstream in(m_inputFileName.c_str());
		if (in.fail())  { cout << "File not found" <<endl; return false; }

		vector<double> *p=NULL;
		
		bool isFirstLine=true;
		while(getline(in, line)  && in.good() )
		{
			vector<string> strv=split(line,",");
			if (!isFirstLine&&m_inputDataColNum!=(strv.size())) {
				return false;
			}
			m_inputDataColNum=strv.size();
			for (int i=0;i<strv.size();i++)
			{
				if ((isFirstLine&&collabel&&rowlabel&&i!=0)||(isFirstLine&&collabel&&!rowlabel))
				{
					m_varColNames.push_back(strv[i]);
				}
				if (!isFirstLine||(isFirstLine&&!collabel)) {
					if (i==0) p=new vector<double>();
					if (i==0&&rowlabel) {
						m_varRowNames.push_back(strv[i]);
					}
					if (i!=0) {
						if (strv[i].empty()) p->push_back(0);
						else{
							double v=str2float(strv[i]);
							p->push_back(v);
						}
					}

				}
			}
			m_inputDataRowNum++;
			int x=rowlabel?m_inputDataColNum-1:m_inputDataColNum;
			if (p!=NULL&&p->size()==x) {
				inputdata.push_back(*p);

			}else if(p!=NULL&&p->size()!=x) return false;

			isFirstLine=false;
		}
		in.close();
	}else return false;
	if (collabel) m_inputDataRowNum--;
	if (rowlabel) m_inputDataColNum--;
	m_inputDataMatrix=new double*[m_inputDataRowNum];

	for (int i=0; i<m_inputDataRowNum; i++) {
		m_inputDataMatrix[i]=new double[m_inputDataColNum];
		for (int j=0; j<m_inputDataColNum; j++) {
			m_inputDataMatrix[i][j]=inputdata[i][j];
		}
	}
	inputdata.clear();
	return true;
}

void MINE::releaseInputMatrix()
{
	if (m_inputDataMatrix!=NULL)
	{
		for (int i=0; i<m_inputDataRowNum; i++) 
			delete []m_inputDataMatrix[i];
	}
	delete []m_inputDataMatrix;
	m_inputDataMatrix=NULL;	
}

bool MINE::parserArgs( int argc, char **argv )
{	
	invocation_name = argv[0];
	const Arg_parser::Option options[] =
	{
		{ 'V', "version",  Arg_parser::no    },
		{ 'a', "alpha",    Arg_parser::yes   },
		{ 'c', "clumps",   Arg_parser::yes },
		{ 'h', "help",     Arg_parser::no    },
		{ 'i', "input",    Arg_parser::yes   },
        { 'l', "inputlabel",    Arg_parser::yes   },
		{ 'o', "output",   Arg_parser::yes   },
        { 'L', "outputlabel",    Arg_parser::yes   },
		{ 'm', "master",    Arg_parser::yes    },
		{ 'A', "allPairs",    Arg_parser::no    },
		{ 'b', "pairsBetween",    Arg_parser::yes    },
		{ 'p', "onePair",    Arg_parser::yes    }
	};

	const Arg_parser parser( argc, argv, options );
	if( parser.error().size() )				// bad option
	{ show_error( parser.error().c_str(), 0, true ); return false; }

	for( int argind = 0; argind < parser.arguments(); ++argind )
	{
		const int code = parser.code( argind );
		string stralpha = parser.argument( argind );
		vector<string> strsplit;
		if( !code ) break;				// no more options
		switch( code )
		{
		case 'V': show_version(); return false;
		case 'i': 
			m_inputFileName=stralpha;
			break;	
		case 'a': 
			param.alpha=str2float(stralpha);
			break;				
		case 'b': 
			m_AnalysisStyle=TwoSets;			
			m_betweenId=str2int(stralpha);
			break;
		case 'c':
			param.c=str2float(stralpha);
			break;				// example, do nothing
		case 'h': show_help(); return false;
		case 'o': 
			m_outputFileName=stralpha;
			break;				// example, do nothing
		case 'm': 
			m_AnalysisStyle=MasterVariable;	
			m_masterId=str2int(stralpha);
			break;
		case 'A': m_AnalysisStyle=AllParis;	 break;
        case 'l':
            m_inputLabel=str2int(stralpha);
            break;
        case 'L':
            m_outputLabel=str2int(stralpha);
            break;
		case 'p':
			strsplit=split(stralpha,",");
			if (strsplit.size()==2)
			{
				m_AnalysisStyle=OneParis;
				m_onepair1=str2int(strsplit[0]);
				m_onepair2=str2int(strsplit[1]);
			}
			break;
		default : internal_error( "uncaught option" );return false; break;
		}
	} // end
	return true;
}

void MINE::printResult(){
    if (m_ResultsArrayLen>0) {
        std::cout<<"var1,var2,mic,mev,mcn,mas"<<endl;
        switch (m_AnalysisStyle) {
            case OneParis:
                std::cout<<m_onepair1<<","<<m_onepair2<<"," <<m_Results[0].mic<<","<<m_Results[0].mev<<","<<m_Results[0].mcn<<","<<m_Results[0].mas<<endl;
                break;
            case AllParis:
                for (int i=0; i<m_inputDataRowNum; i++) {
                    for (int j=i+1; j<m_inputDataRowNum; j++) {
                        int id=i*(m_inputDataRowNum-1)-i*(i-1)/2+j-i-1;
                        std::cout<<i<<","<<j<<"," <<m_Results[id].mic<<","<<m_Results[id].mev<<","<<m_Results[id].mcn<<","<<m_Results[id].mas<<endl;
                    }
               }
                break;
            case MasterVariable:
                for (int i=0; i<m_inputDataRowNum; i++) {
                    if (i==m_masterId) {
                        continue;
                    }else{
                        std::cout<<m_masterId<<","<<i<<"," <<m_Results[i].mic<<","<<m_Results[i].mev<<","<<m_Results[i].mcn<<","<<m_Results[i].mas<<endl;
                    }
                }
                break;
            case TwoSets:
                for (int i=0; i<m_betweenId; i++) {
                    for (int j=m_betweenId; j<m_inputDataRowNum; j++) {
                        int id=i*m_betweenId+j-m_betweenId;
                        std::cout<<i<<","<<j<<"," <<m_Results[id].mic<<","<<m_Results[id].mev<<","<<m_Results[id].mcn<<","<<m_Results[id].mas<<endl;
                    }
                }
            default:
                break;
        }
    }
}

bool MINE::exportResult()
{
    cout << "begin writing result file;"<<endl;
    if (m_ResultsArrayLen>0) {
        if (m_varRowNames.size()>0&&m_outputLabel==1) {
            ofstream outputfile(m_outputFileName.c_str());
            if (outputfile.is_open()) {
                outputfile<<"var1,var2,mic,mev,mcn,mas"<<endl;
                switch (m_AnalysisStyle) {
                    case OneParis:
                        outputfile<<m_varRowNames[m_onepair1]<<","<<m_varRowNames[m_onepair2]<<"," <<m_Results[0].mic<<","<<m_Results[0].mev<<","<<m_Results[0].mcn<<","<<m_Results[0].mas<<endl;
                        break;
                    case AllParis:
                        for (int i=0; i<m_inputDataRowNum; i++) {
                            for (int j=i+1; j<m_inputDataRowNum; j++) {
                                int id=i*(m_inputDataRowNum-1)-i*(i-1)/2+j-i-1;
                                outputfile<<m_varRowNames[i]<<","<<m_varRowNames[j]<<"," <<m_Results[id].mic<<","<<m_Results[id].mev<<","<<m_Results[id].mcn<<","<<m_Results[id].mas<<endl;
                            }
                        }
                        break;
                    case MasterVariable:
                        for (int i=0; i<m_inputDataRowNum; i++) {
                            if (i==m_masterId) {
                                continue;
                            }else{
                                outputfile<<m_varRowNames[m_masterId]<<","<<m_varRowNames[i]<<"," <<m_Results[i].mic<<","<<m_Results[i].mev<<","<<m_Results[i].mcn<<","<<m_Results[i].mas<<endl;
                            }
                        }
                        break;
                    case TwoSets:
                        for (int i=0; i<m_betweenId; i++) {
                            for (int j=m_betweenId; j<m_inputDataRowNum; j++) {
                                int id=i*m_betweenId+j-m_betweenId;
                                outputfile<<m_varRowNames[i]<<","<<m_varRowNames[j]<<"," <<m_Results[id].mic<<","<<m_Results[id].mev<<","<<m_Results[id].mcn<<","<<m_Results[id].mas<<endl;
                            }
                        }
                    default:
                        break;
                }
            }else {cout<< "Output file can not be created!"; return false;}
            outputfile.close();
            
        }else{
            ofstream outputfile(m_outputFileName.c_str());
            if (outputfile.is_open()) {
                outputfile<<"var1,var2,mic,mev,mcn,mas"<<endl;
                switch (m_AnalysisStyle) {
                    case OneParis:
                        outputfile<<m_onepair1<<","<<m_onepair2<<"," <<m_Results[0].mic<<","<<m_Results[0].mev<<","<<m_Results[0].mcn<<","<<m_Results[0].mas<<endl;
                        break;
                    case AllParis:
                        for (int i=0; i<m_inputDataRowNum; i++) {
                            for (int j=i+1; j<m_inputDataRowNum; j++) {
                                int id=i*(m_inputDataRowNum-1)-i*(i-1)/2+j-i-1;
                                outputfile<<i<<","<<j<<"," <<m_Results[id].mic<<","<<m_Results[id].mev<<","<<m_Results[id].mcn<<","<<m_Results[id].mas<<endl;
                            }
                        }
                        break;
                    case MasterVariable:
                        for (int i=0; i<m_inputDataRowNum; i++) {
                            if (i==m_masterId) {
                                continue;
                            }else{
                                outputfile<<m_masterId<<","<<i<<"," <<m_Results[i].mic<<","<<m_Results[i].mev<<","<<m_Results[i].mcn<<","<<m_Results[i].mas<<endl;
                            }
                        }
                        break;
                    case TwoSets:
                        for (int i=0; i<m_betweenId; i++) {
                            for (int j=m_betweenId; j<m_inputDataRowNum; j++) {
                                int id=i*m_betweenId+j-m_betweenId;
                                outputfile<<i<<","<<j<<"," <<m_Results[id].mic<<","<<m_Results[id].mev<<","<<m_Results[id].mcn<<","<<m_Results[id].mas<<endl;
                            }
                        }
                    default:
                        break;
                }
            }else {cout<< "Output file can not be created!"; return false;}
            outputfile.close();
            
        }
    }else return false;
    cout << "end all!"<<endl;
    return  true;
    
    
}