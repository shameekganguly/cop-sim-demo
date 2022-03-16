// event logger
#ifndef EVENT_LOGGER_H
#define EVENT_LOGGER_H

#include <iostream>
#include <string>

/*  THIS IS WIP CODE */

enum EventLogLevel {
	LOG_DEBUG = 0, // to be used for debugging, probably in debug mode. logs everything.
	LOG_DIAGNOSE, // to be used for logging diagnostic data to figure out algorithmic malfunctions.
	LOG_PROD // to be used in production. Only logs exceptions.
};

struct LogEvent {
	std::string mssg;
	EventLogLevel level;
	system_clock::time_point logtime;
	LogEvent(const std::string inmssg, EventLoglevel inlevel)
	mssg(inmssg), level(inlevel)
	{
		logtime = system_clock::now();
	}
};

// threaded event logger
class EventLogger {
public:
	// constructor
	EventLogger (const std::string& fname, const EventLogLevel level, unsigned long intervalms = 500)
	: _allowed_level(level) {
		_f_log_cout = false;
		if (fname.empty()) {
			_f_log_cout = true;
		} else {
			_log_file.open(fname);
			if(!_log_file.is_open()) {
				_f_log_cout = true;
			}
		}

		const unsigned long MIN_LOG_FREQ = 100;
		if(intervalms < MIN_LOG_FREQ) {
			std::cerr << "Event Logger tried to start with very high frequency: " << intervalms <<;
			std::cerr << "Resetting to " << MIN_LOG_FREQ << "ms instead" << std::endl;
			intervalms = MIN_LOG_FREQ;
		}
		_log_interval = intervalms;

		resizeQueue();
	}

	// set log level
	void allowedLogLevelIs (const EventLogLevel level) {
		if(_allowed_level == level) return;
		_allowed_level = level;
		resizeQueue();
	}

	// clear message queue
	void clearMessageQueue () {
		uint num_messages_to_clear = _events.size();

	}

	// log filename
	std::string _log_fname;

	// start time
	system_clock::time_point _t_start;

	// interval how often the queue will be processed for file/ stream write. in milliseconds
	unsigned long _log_interval;

	// log level
	EventLogLevel _allowed_level;

private:
	void resizeQueue() {
		// set queue size based on current log level, interval, and expected log frequency
		unsigned long message_queue_size;
		const unsigned long MIN_QUEUE_SIZE = 100;
		switch(_allowed_level) {
			case LOG_DEBUG:
				// expect upto 20 messages per millisecond
				message_queue_size = 20 * _log_interval;
				break;
			case LOG_DIAGNOSE:
				// expect upto 5 messages per millisecond
				message_queue_size = 5 * _log_interval;
				break;
		    case LOG_PROD:
		    	// expect upto 10 messages per second (0.01)
		    	message_queue_size = max(MIN_QUEUE_SIZE, (10*_log_interval)/1000);
		    	break;
		}
		if(message_queue_size > _events.size()) {
			_events.reserve()
		}
	}

	// size of the log queue. all messages are discarded when the log is full.
	uint _queue_size;

	// log queue
	std::queue<LogEvent> _events;

	// should log to cout instead?
	bool _f_log_cout;

	// log file
	std::fstream _log_file;

	// thread
	std::thread _log_thread;
};

#endif // EVENT_LOGGER_H