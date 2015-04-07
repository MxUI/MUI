/*
 * lib_uri.h
 *
 *  Created on: Mar 14, 2014
 *      Author: ytang
 */

#ifndef LIB_URI_H_
#define LIB_URI_H_

#include "util.h"

namespace mui {

class uri {
public:
	uri(const std::string& url_s) {
		parse(url_s);
	}
	uri(const char url_c[]) {
		parse(url_c);
	}
	const std::string& protocol() const { return protocol_; }
	const std::string& host()     const { return host_; }
	const std::string& path()     const { return path_; }

	uri( const uri &another ) = delete;
	uri& operator = ( const uri &another ) = delete;
private:
    void parse(const std::string& url_s) {
    	// "__protocol__://__host__/__path__"
	std::size_t prot_end = url_s.find("://");
	protocol_ = url_s.substr(0,prot_end);
	std::size_t host_end = url_s.find("/",prot_end+3);
	host_ = url_s.substr(prot_end+3,host_end-prot_end-3);
	path_ = url_s.substr(host_end+1);

	std::transform(protocol_.begin(), protocol_.end(), protocol_.begin(), ::tolower);
	std::transform(host_.begin(), host_.end(), host_.begin(), ::tolower);
    }

    std::string protocol_, host_, path_;
};

}

#endif /* LIB_URI_H_ */
