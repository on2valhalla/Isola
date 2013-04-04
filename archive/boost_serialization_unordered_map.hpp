//
//  boost_serialization_unordered_map.h
//  isola
//
//  Created by Jason Mann on 3/31/13.
//  Copyright (c) 2013 Jason Mann. All rights reserved.
//

#ifndef isola_boost_serialization_unordered_map_h
#define isola_boost_serialization_unordered_map_h

#include <unordered_map>

#include <boost/config.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/collections_load_imp.hpp>
#include <boost/serialization/split_free.hpp>

namespace boost {
	namespace serialization {
		
		template<class Archive, class Type, class Key, class Hash, class
		Compare, class Allocator >
		inline void save(
						 Archive & ar,
						 const std::unordered_map<Key, Type, Hash, Compare, Allocator> &t,
						 const unsigned int /* file_version */
						 ){
			boost::serialization::stl::save_collection<
			Archive,
			std::unordered_map<Key, Type, Hash, Compare, Allocator>
			>(ar, t);
		}
		
		template<class Archive, class Type, class Key, class Hash, class
		Compare, class Allocator >
		inline void load(
						 Archive & ar,
						 std::unordered_map<Key, Type, Hash, Compare, Allocator> &t,
						 const unsigned int /* file_version */
						 ){
			boost::serialization::stl::load_collection<
			Archive,
			std::unordered_map<Key, Type, Hash, Compare, Allocator>,
			boost::serialization::stl::archive_input_seq<
			Archive, std::unordered_map<Key, Type, Hash, Compare,
			Allocator> >,
			
			boost::serialization::stl::no_reserve_imp<std::unordered_map<
			Key, Type, Hash, Compare, Allocator
			>
			>
			>(ar, t);
		}
		
		// split non-intrusive serialization function member into separate
		// non intrusive save/load member functions
		template<class Archive, class Type, class Key, class Hash, class
		Compare, class Allocator >
		inline void serialize(
							  Archive & ar,
							  std::unordered_map<Key, Type, Hash, Compare, Allocator> &t,
							  const unsigned int file_version
							  ){
			boost::serialization::split_free(ar, t, file_version);
		}
		
	} // serialization
} // namespace boost


#endif
