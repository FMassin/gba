#define SEISCOMP_COMPONENT TEST

#include <seiscomp3/logging/log.h>
#include <seiscomp3/client/streamapplication.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/io/records/mseedrecord.h>
#include <seiscomp3/processing/eewamps/processor.h>
#include <string>
#include <gba.h>


using namespace std;
using namespace Seiscomp;


class SCGbA : public Client::StreamApplication {
	public:
		SCGbA(int argc, char** argv) : Client::StreamApplication(argc, argv) {
			setMessagingEnabled(false);
			setDatabaseEnabled(true, true);
			setLoadStationsEnabled(true);
			setLoggingToStdErr(true);
			setRecordDatatype(Array::FLOAT);
		}


		void createCommandLineDescription() {
			Client::Application::createCommandLineDescription();

			commandline().addGroup("Streams");
			commandline().addOption("Streams", "streams-allow,a", "Stream IDs to allow separated by comma", &_allowString);
			commandline().addOption("Streams", "streams-deny,r", "Stream IDs to deny separated by comma", &_denyString);
			commandline().addOption("Streams", "dump", "Dump all processed streams as mseed to stdout");
		}


		bool validateParameters() {
			if ( !StreamApplication::validateParameters() )
				return false;

			if ( !_allowString.empty() ) {
				std::vector<std::string> tokens;
				Core::split(tokens, _allowString.c_str(), ",");
				for ( size_t i = 0; i < tokens.size(); ++i )
					_eewProc.addAllowRule(tokens[i]);
			}

			if ( !_denyString.empty() ) {
				std::vector<std::string> tokens;
				Core::split(tokens, _denyString.c_str(), ",");
				for ( size_t i = 0; i < tokens.size(); ++i )
					_eewProc.addDenyRule(tokens[i]);
			}
			if ( !isInventoryDatabaseEnabled() )
				setDatabaseEnabled(false, false);

			return true;
		}


		bool init() {
			_first = true;
			if ( !StreamApplication::init() )
				return false;

			//Create test picks
			_pickabk = new DataModel::Pick("TESTPICK_ABK");
			Core::Time ptime(1996,8,10,18,12,21,840000);
			Core::Time now = Core::Time::GMT();
			DataModel::CreationInfo ci;
			ci.setCreationTime(now);
			ci.setAgencyID("GbA");
			_pickabk->setCreationInfo(ci);
			_pickabk->setTime(DataModel::TimeQuantity(ptime));
			_pickabk->setWaveformID(DataModel::WaveformStreamID("BO","ABK","","HGZ","SomeURI/ABK"));
			_pickabk->setPhaseHint(DataModel::Phase("P"));

			_pickabo = new DataModel::Pick("TESTPICK_ABO");
			ci.setCreationTime(now);
			ci.setAgencyID("GbA");
			_pickabo->setCreationInfo(ci);
			_pickabo->setTime(DataModel::TimeQuantity(ptime));
			_pickabo->setWaveformID(DataModel::WaveformStreamID("BO","ABO","","HGZ","SomeURI/ABO"));
			_pickabo->setPhaseHint(DataModel::Phase("P"));


			Processing::EEWAmps::Config eewCfg;
			eewCfg.dumpRecords = commandline().hasOption("dump");

			eewCfg.baseLineCorrectionBufferLength = 60.0;
			eewCfg.taperLength = 10.0;
			eewCfg.gba.enable = true;
			eewCfg.gba.bufferSize = Core::TimeSpan(30,0);
			eewCfg.gba.cutOffTime = Core::TimeSpan(20,0);

			// Convert to velocity only
			eewCfg.wantSignal[Processing::WaveformProcessor::MeterPerSecond] = true;

			_eewProc.setConfiguration(eewCfg);
			_eewProc.setGbACallback(boost::bind(&SCGbA::handleFilterBank, this, _1, _2, _3, _4, _5, _6));
			_eewProc.setInventory(Client::Inventory::Instance()->inventory());

			if ( !_eewProc.init(configuration()) )
				return false;

			_eewProc.showConfig();
			_eewProc.showRules();
			_eewProc.subscribeToChannels(recordStream(), Core::Time::GMT());

			// Initialise magnitude and distance samples; pdf array; marginal
			// arrays
			_nms = 31;
			_magmin = 2.0;
			_dmag = 0.2;
			_nrs = 21;
			_rmin = 0.0;
			_dr = 0.1;
			_msamples = new double[_nms];
			for(int i=0;i<_nms;i++)
				_msamples[i] = _magmin + i*_dmag;
			_rsamples = new double[_nrs];
			for(int i=0;i<_nrs;i++)
				_rsamples[i] = _rmin + i*_dr;

			_mmarg = new double[_nms];
			_rmarg = new double[_nrs];
			_pdf = new double[_nms*_nrs];

			return true;
		}

		void handleRecord(Record *rec) {
			RecordPtr tmp(rec);
			_eewProc.feed(rec);
			Core::Time start(rec->startTime());
			if ((double)start > (double)_pickabk->time().value() && _first){
				_eewProc.feed(_pickabk);
				_eewProc.feed(_pickabo);
				SEISCOMP_DEBUG("Sending pick at %s",start.iso().c_str());
				_first = false;
			}
		}

		void get_moments(double *pdf, double *samples, int ndata, double ds,
				         double &mean, double &max, double &sigma){
			int i, imax=0;
			double pdfmax = std::numeric_limits<double>::min();
			double tmp1, tmp2;

			// normalize pdf
			double evd = 0;
			for(i=0;i<ndata-1;i++)
				evd = evd + (pdf[i] + pdf[i+1])/2.*ds;
			for(i=0;i<ndata;i++)
				pdf[i] = pdf[i] / evd;

			mean = 0;
			for(i=0; i<ndata-1; i++){
				mean += (pdf[i]*samples[i] + pdf[i+1]*samples[i+1])/2.*ds;
				if (pdf[i] > pdfmax){
					pdfmax = pdf[i];
					imax = i;
				}
			}
			if(pdf[ndata-1] > pdfmax)
				imax = ndata-1;
			max = samples[imax];

			sigma=0;
			for(i=0; i<ndata; i++){
				tmp1 = (pdf[i]*(samples[i]-mean)*(samples[i]-mean))/2.;
				tmp2 = (pdf[i+1]*(samples[i+1]-mean)*(samples[i+1]-mean))/2.;
				sigma += (tmp1+tmp2)*ds;
			}
		}

		void evaluate(){
			if (_mags.empty())
				return;

			double mhat, rhat, pdfmax, mbar, rbar, m_sig, r_sig;
			int i, imhat, irhat;
			bool first = true;
			double _mmarg_tmp[_nms], _rmarg_tmp[_nrs];

			memset(_mmarg,0,_nms*sizeof(double));
			memset(_rmarg,0,_nrs*sizeof(double));
			for(Estimates::iterator it = _mags.begin(); it != _mags.end(); ++it){
				GbA* est = it->second;
				memset(_mmarg_tmp,0,_nms*sizeof(double));
				memset(_rmarg_tmp,0,_nrs*sizeof(double));
				memset(_pdf,0,_nms*_nrs*sizeof(double));
				est->get_pdf(_pdf,_nms,_nrs,&_mmarg_tmp[0],_nms,
						     &_rmarg_tmp[0],_nrs);
				for(i=0; i<_nms;i++){
					if(first){
						_mmarg[i] = _mmarg_tmp[i];
					}else{
						_mmarg[i] *= _mmarg_tmp[i];
					}
				}
				for(i=0; i<_nrs;i++){
					if(first){
						_rmarg[i] = _rmarg_tmp[i];
					}else{
						_rmarg[i] *= _rmarg_tmp[i];
					}
				}
				first=false;
			}
			get_moments(_mmarg,_msamples,_nms,_dmag,mbar,mhat,m_sig);
			get_moments(_rmarg,_rsamples,_nrs,_dr,rbar,rhat,r_sig);
			cout << _latest_data.toString("%FT%T.%fZ");
			cout << "; " << Core::Time::GMT().toString("%FT%T.%fZ");
			cout << "; Mhat: " << mhat << "; Mbar: " << mbar <<"; Sigma^2: " << m_sig;
			cout << "; Rhat: " << rhat << "; Rbar: " << rbar <<"; Sigma^2: " << r_sig << endl;
			return;
		}

		void handleFilterBank(const Processing::EEWAmps::BaseProcessor *proc,
		                      std::string pickID, double *amps,
		                      const Core::Time &max, const Core::Time &ptime,
		                      const Core::Time &maxevaltime){
			double mean[2], cov[2][2];
			DataModel::WaveformStreamID wid = proc->waveformID();
			std::string stream(wid.networkCode()+"."+wid.stationCode()+"."+wid.locationCode());
			Estimates::iterator it = _mags.find(stream);
			if(it == _mags.end()){
				SEISCOMP_DEBUG("Init GbA processing for stream %s",stream.c_str());
				_mags[stream] = new GbA(_eewProc.configuration().gba.passbands.size(),30);
				_mags[stream]->init("/home/behry/workspace/gutenberg_algorithm/data/GbA_training.nc");
			}
			_mags[stream]->process(amps,_eewProc.configuration().gba.passbands.size(),
					              (double)max - (double)ptime,proc->usedComponent());
			//_mags[stream]->get_mean_cov(mean,cov);
			if((double)max > (double)_latest_data)
				_latest_data = (double)max;

			/*cout << proc->usedComponent() << "; " << proc->streamID() <<"; ";
			cout << ptime.iso() << "; " << max.iso() << "; " << (double)max - (double)ptime << "; ";
			cout << mean[0] << ";" << mean[1];
			cout << endl;*/
			evaluate();
		}


	private:
		typedef map<std::string, GbA*> Estimates;
		std::string                        _allowString, _denyString;
		Processing::EEWAmps::Processor     _eewProc;
		DataModel::Pick                    *_pickabk, *_pickabo;
		Estimates                          _mags;
		int                                _nms, _nrs;
		double                             _magmin, _dmag;
		double                             _rmin, _dr;
		double                             *_mmarg, *_rmarg;
		double                             *_pdf;
		double                             *_msamples, *_rsamples;
		bool                               _first;
		Core::Time                         _latest_data;
};


int main(int argc, char **argv) {
	return SCGbA(argc, argv)();
}
