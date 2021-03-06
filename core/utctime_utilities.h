/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <iosfwd>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <cmath>
#include "core_serialization.h"

namespace shyft {
	namespace core {
        /** \brief utctime
        * basic types for time handling
        * we use linear time, i.e. time is just a number
        * on the timeaxis, utc. Timeaxis zero is at 1970-01-01 00:00:00 (unix time).
        * resolution is 1 utctimespan, currently micro seconds.
        * also define min max and no_utctime
        *
        * The advantage with the definition is that it is well defined and commonly known in all platforms
        */
        using utc_clock = std::chrono::system_clock;
        using utctimespan = std::chrono::microseconds;
        using seconds = std::chrono::seconds;
        inline utctimespan deltahours(int h) { return std::chrono::duration_cast<utctimespan>(std::chrono::hours(h)); }
        inline utctimespan deltaminutes(int m) { return std::chrono::duration_cast<utctimespan>(std::chrono::minutes(m)); }

        using utctime = utctimespan;//std::chrono::time_point<utc_clock,utctimespan>;// clock::time_point<chrono::milliseconds>;
        constexpr const utctime utctime_0 = utctime{seconds{0}};// epoch, 1970.01.01 utc
        constexpr const utctime max_utctime	= utctime::max();	/// max 64bit int
        constexpr const utctime min_utctime	= utctime::min()+utctimespan{1};
        constexpr const utctime no_utctime = utctime::min();

        /** \brief current utctime
        *  \return current systemclock utctime
        */
        inline utctime utctime_now() {return std::chrono::time_point_cast<utctimespan>(utc_clock::now()).time_since_epoch(); }

        inline bool is_valid(utctime t) {return t != no_utctime;}

        /** \brief computes floor of t vs utctimespan dt
        *
        * If dt is 0, t is returned
        * if dt< 0 then it computes ceil, (hmm)
        * \return floor of t vs dt as explained above
        */
        inline utctime floor(utctime t, utctimespan dt) noexcept {
            int64_t den = dt.count();
            if (den == 0)
                return t;
            int64_t num = t.count();
            if (0 < (num^den))
                return utctime(utctimespan(den*(num / den)));
            auto r = lldiv(num, den);
            return utctime(utctimespan(r.rem ? den*(r.quot - 1) : den*r.quot));
        }

        utctime utctime_floor(utctime t, utctimespan dt) noexcept; // expose to python

        utctime create_from_iso8601_string(const std::string&s);

        inline double to_seconds(const utctimespan &dt) noexcept { return double(dt.count()) / utctimespan{seconds{1}}.count(); }
        inline int64_t to_seconds64(const utctimespan&dt ) noexcept {return dt.count()/utctimespan{seconds{1}}.count();}
        inline utctimespan from_seconds(double sec) noexcept { return utctimespan{ int64_t(round(utctimespan::period::den*sec / utctimespan::period::num)) }; } // NAN,INF?
        constexpr inline utctime utctime_from_seconds64(int64_t x) noexcept {return utctime{seconds(x)};}


        struct calendar;

        enum class trim_policy:int8_t {TRIM_IN, TRIM_OUT};

        /** \brief utcperiod is defined
         *  as period on the utctime space, like
         * [start..end>, where end >=start to be valid
         *
         */
		struct utcperiod {
            utcperiod(utctime start, utctime end) noexcept: start(start), end(end) {}
            utcperiod(int64_t start, int64_t end) noexcept : start{ seconds(start) }, end{ seconds(end) } {}
            utcperiod() noexcept: start(no_utctime), end(no_utctime) {}
            utctimespan timespan() const noexcept {	return (end - start); }
			bool valid() const noexcept { return start != no_utctime && end != no_utctime && start <= end; }
			bool operator==(const utcperiod &p) const noexcept { return start == p.start &&  end == p.end; }
			bool operator!=(const utcperiod &p) const noexcept {return ! (*this==p);}
			bool contains(utctime t) const noexcept {return is_valid(t)&&valid()?t>=start && t<end:false;}
            bool contains(int64_t t) const noexcept {return contains(seconds{t});}
			bool contains(const utcperiod& p) const noexcept {return valid()&&p.valid()&& p.start>=start &&p.end <=end;}
			bool overlaps(const utcperiod& p) const noexcept {return ( (p.start >= end) || (p.end <= start) )?false:true; }

            utcperiod trim(const calendar &c, utctimespan deltaT, trim_policy tp) const;
            utcperiod trim(const calendar &c, int64_t deltaT, trim_policy tp) const {return trim(c,seconds{deltaT},tp);}
            int64_t diff_units(const calendar &c, utctimespan deltaT) const;
            int64_t diff_units(const calendar &c, int64_t deltaT) const {return diff_units(c,seconds{deltaT});}
            int64_t diff_units(const calendar &c, utctimespan deltaT, utctimespan &remainder) const;
            int64_t diff_units(const calendar &c, int64_t deltaT, utctimespan &remainder) const {return diff_units(c,seconds{deltaT},remainder);}

			utctime start;
			utctime end;
			std::string to_string() const;
            friend std::ostream& operator<<(std::ostream& os, const utcperiod& p);
            x_serialize_decl();
		};
		inline bool is_valid(const utcperiod &p) noexcept {return p.valid();}
		/** returns intersection of a and b
         *  ensures that if there is an a.overlaps(b), the intersection is returned, where .timespan() >0
         *  otherwise an empty not .valid() period is returned.
         *  \note: that intersection of 0 is not defined as intersection
         */
        inline utcperiod intersection(const utcperiod&a, const utcperiod& b) noexcept {
			if (a.overlaps(b)) {
				return utcperiod(std::max(a.start, b.start), std::min(a.end, b.end));
			}
			return utcperiod();
        };

        namespace time_zone {
            using namespace std;

            /**\brief time_zone handling, basically just a wrapper that hides
            * the fact we are using boost::date_time for getting/providing tz info
            */
            template<typename tz>
            struct tz_info {
                typedef tz tz_type_t;
                string name() const {return "UTC";}
                utctimespan base_offset() const {return utctimespan(0);}
                utctimespan utc_offset(utctime /*t*/) const {return utctimespan(0);}
                bool is_dst(utctime t) const {return false;}
                utctime dst_start(int /* y*/) const {return no_utctime;}
                utctime dst_end(int /*y*/) const {return no_utctime;}
            };


            /**\brief The tz_table is a table driven approach where each year in a specified range have
             * a dst information containing start/end, and the applied dst offset.
             * This approach allows to have historical correct dst-rules, at minor space/time overhead
             * E.g. Norway, summertime rules are changed to EU defs. in 1996
             * first time applied was 1916, then 1943-1945, then 1959-1965, then 1980-to 1996,
             * where eu rules was introduced.
             */
            struct tz_table {
                int64_t start_year;
                string tz_name;
                vector<utcperiod> dst;
                vector<utctimespan> dt;

                /**\brief construct a tz_table using a information from provided Tz adapter
                 *\tparam Tz a type that provids dst_start(year),dst_end(year), dst_offset(year)
                 *\param tz const ref to provide tz information that will be used to construct a table driven interpretation
                 *\param start_year default 1905 (limit of posix time is 1901) for which the dst table is constructed
                 *\param n_years default 200 giving range from 1905 to 2105 with dst info.
                 */
                template<typename Tz>
                tz_table(const Tz& tz ,int start_year=1905,size_t n_years=200):start_year(start_year) {
                    for(int y=start_year;y<int(start_year+n_years);++y) {
                        dst.emplace_back(tz.dst_start(y),tz.dst_end(y));
                        dt.push_back(tz.dst_offset(y));
                    }
                    tz_name=tz.name();
                }
                /**\brief construct a simple dst infotable with no dst, just tz-offset
                * suitable for non-dst time-zones and data-exchange.
                * \param dt of type utctimespan, positive for tz east of GMT
                */
                explicit tz_table(utctimespan dt):start_year(0) {
                    if(dt != utctimespan{0}) {
                        char s[100];sprintf(s,"UTC%+02d",int(dt/deltahours(1)));
                        tz_name=s;
                    } else {
                        tz_name="UTC";
                    }
                }
                tz_table():start_year(0),tz_name("UTC"){}
                inline bool is_dst() const {return dst.size()>0;}
                string name() const {return tz_name;}
                utctime dst_start(int year) const {return is_dst()?dst[year-start_year].start:no_utctime;}
                utctime dst_end (int year) const {return is_dst()?dst[year-start_year].end:no_utctime;}
                utctimespan dst_offset(utctime t) const ;
                x_serialize_decl();
            };

            /**\brief a table driven tz_info, using the tz_table implementation */
            template<>
            struct tz_info<tz_table> {
                utctimespan base_tz;
                tz_table tz;
                tz_info():base_tz(0) {}// serialization
                tz_info(utctimespan base_tz):base_tz(base_tz),tz(base_tz) {}
                tz_info(utctimespan base_tz,const tz_table&tz):base_tz(base_tz),tz(tz) {}
                string name() const {return tz.name();}
                utctimespan base_offset() const {return base_tz;}
                utctimespan utc_offset(utctime t) const {return base_tz + tz.dst_offset(t);}
                bool is_dst(utctime t) const {return tz.dst_offset(t)!=utctimespan(0);}
                bool is_any_dst() const {return tz.is_dst();}
                utctime dst_start(int year) const {return tz.dst_start(year);}
                utctime dst_end(int year) const {return tz.dst_end(year);}
                x_serialize_decl();
            };

            typedef tz_info<tz_table> tz_info_t;///< tz_info type most general, supports all kinds of tz, at a minor extra cost.
            typedef shared_ptr<tz_info_t> tz_info_t_;///< convinience, the shared ptr. version

            /** \brief time zone database class that provides shared_ptr of predefined tz_info_t objects */
            struct tz_info_database {

                /** \brief load from compile time iso db as per boost 1.60 */
                void load_from_iso_db();

                /** \brief load from file that contains all descriptions, ref. boost::date_time for format */
                void load_from_file(const string& filename);

                /** \brief add one entry, using a specified region_name like Europe/Copenhagen, and a posix description string, ref boost::date_time for spec */
                void add_tz_info(string region_name,string posix_tz_string);

                /** \brief returns a shared_ptr to tz_info_t given time-zone region name like Europe/Copenhagen */
                shared_ptr<tz_info_t> tz_info_from_region(const string &region_name) const {
                    auto f=region_tz_map.find(region_name);
                    if( f!= region_tz_map.end()) return f->second;
                    throw runtime_error(string("tz region '")+region_name + string("' not found"));
                }

                /** \brief returns a shared_ptr to tz_info_t given time-zone name like CET */
                shared_ptr<tz_info_t> tz_info_from_name(const string &name) const {
                    auto f=name_tz_map.find(name);
                    if( f!= name_tz_map.end()) return f->second;
                    throw runtime_error(string("tz name '")+name + string("' not found"));
                }
                vector<string> get_region_list() const {
                    vector<string> r;r.reserve(region_tz_map.size());
                    for(const auto& c:region_tz_map)
                        r.push_back(c.first);
                    return r;
                }
                vector<string> get_name_list() const {
                    vector<string> r;r.reserve(region_tz_map.size());
                    for(const auto& c:name_tz_map)
                        r.push_back(c.first);
                    return r;
                }

                map<string,shared_ptr<tz_info_t>> region_tz_map;///< map from Europe/Copenhagen to tz
                map<string,shared_ptr<tz_info_t>> name_tz_map;///< map from CET to tz (same tz as Europe/Copenhagen)
            };

        }
        /** \brief YMDhms, simple structure that contains calendar coordinates.
         * Contains year,month,day,hour,minute, second,
         * for ease of constructing utctime.
         *\note the constructor do a range check for Y M D h m s, and throws if fail.
         *
         */
		struct YMDhms {
            static const int YEAR_MAX= 9999;
            static const int YEAR_MIN=-9999;
			YMDhms():year(0), month(0), day(0), hour(0), minute(0), second(0),micro_second(0) {}
			YMDhms(int Y, int M=1, int D=1, int h=0, int m=0, int s=0,int us=0) : year(Y), month(M), day(D), hour(h), minute(m), second(s),micro_second(us)  {
                if(!is_valid())
                    throw std::runtime_error("calendar coordinates failed simple range check for one or more item:"+std::to_string(Y));
			}

			int year; int month; int day; int hour; int minute; int second;int micro_second;
			///< just check that YMDhms are within reasonable ranges,\note it might still be an 'invalid' date!
			bool is_valid_coordinates() const {return !(year<YEAR_MIN || year>YEAR_MAX || month<1 || month>12 ||day<1 || day>31 ||hour<0 || hour>23 || minute<0 ||minute>59||second<0 ||second>59||micro_second<0 || micro_second>=1000000);}
            ///< if a 'null' or valid_coordinates
			bool is_valid() const { return is_null() || is_valid_coordinates(); }
			bool is_null() const { return year == 0 && month == 0 && day == 0 && hour == 0 && minute == 0 && second == 0 && micro_second==0; }
			bool operator==(YMDhms const& x) const {
                return x.year == year && x.month == month && x.day == day && x.hour == hour
                       && x.minute == minute && x.second == second && x.micro_second==micro_second;
            }
            bool operator!=(YMDhms const&o) const { return !operator==(o); }
			static YMDhms max() {return YMDhms(YEAR_MAX,12,31,23,59,59);}
			static YMDhms min() {return YMDhms(YEAR_MIN,1,1,0,0,0);}
		};
        struct YWdhms {
            int iso_year;
            int iso_week;
            int week_day;
            int hour;
            int minute;
            int second;
            int micro_second;
            YWdhms() :iso_year(0), iso_week(0), week_day(0), hour(0), minute(0), second(0),micro_second(0) {}
            YWdhms(int iso_year,
                int iso_week=1,
                int week_day=1,
                int hour=0,
                int minute=0,
                int second=0,
                int micro_second=0
            ):iso_year(iso_year), iso_week(iso_week), week_day(week_day), hour(hour), minute(minute), second(second),micro_second(micro_second) {
                if (!is_valid())
                    throw std::runtime_error("calendar iso week coordinates failed simple range check for one or more item");
            }
            bool operator==(YWdhms const& x) const {
                return x.iso_year == iso_year && x.iso_week == iso_week && x.week_day == week_day && x.hour == hour
                    && x.minute == minute && x.second == second && x.micro_second==micro_second;
            }
            bool operator!=(YWdhms const&x) const { return !operator==(x); }
            bool is_null() const { return iso_year == 0 && iso_week == 0 && week_day == 0 && hour == 0 && minute == 0 && second == 0 && micro_second==0; }
            ///< just check that YMDhms are within reasonable ranges,\note it might still be an 'invalid' date!
            bool is_valid_coordinates() const { return !(iso_year<YMDhms::YEAR_MIN || iso_year>YMDhms::YEAR_MAX || iso_week < 1 || iso_week>53 || week_day < 1 || week_day>7 || hour < 0 || hour>23 || minute < 0 || minute>59 || second < 0 || second>59 || micro_second<0 || micro_second >=1000000); }
            ///< if a 'null' or valid_coordinates
            bool is_valid() const { return is_null() || is_valid_coordinates(); }
            static YWdhms max() { return YWdhms(YMDhms::YEAR_MAX, 52, 6, 23, 59, 59); }
            static YWdhms min() { return YWdhms(YMDhms::YEAR_MIN, 1, 1, 0, 0, 0); }

        };
        /** \brief Calendar deals with the concept of human calendar.
         *
         * Please notice that although the calendar concept is complete,
         * we only implement features as needed in the core and interfaces.
         *
         * including:
         * -# Conversion between the calendar coordinates YMDhms and utctime, taking  any timezone and DST into account
         * -# Calendar constants, utctimespan like values for Year,Month,Week,Day,Hour,Minute,Second
         * -# Calendar arithmetic, like adding calendar units, e.g. day,month,year etc.
         * -# Calendar arithmetic, like trim/truncate a utctime down to nearest timespan/calendar unit. eg. day
         * -# Calendar arithmetic, like calculate difference in calendar units(e.g days) between two utctime points
         * -# Calendar Timezone and DST handling
         * -# Converting utctime to string and vice-versa
         */
		struct calendar {
			// these do have calendar sematics(could/should be separate typed/enum instad)
			static const utctimespan YEAR;//=seconds{365*24*3600L};
            static const utctimespan QUARTER;// = seconds{3*30*24*3600L};
			static const utctimespan MONTH;// = seconds{30*24*3600L};
			static const utctimespan WEEK;// = seconds{7*24*3600L};
			static const utctimespan DAY;// =  seconds{1*24*3600L};
			static const utctimespan HOUR_3;// = seconds{3 * 3600L};
			// these are just timespan constants with no calendar semantics
			static const utctimespan HOUR;// = seconds{3600L};
			static const utctimespan MINUTE;// = seconds{60L};
			static const utctimespan SECOND;// = seconds{1L};

			static constexpr const int64_t UnixDay = 2440588;///< Calc::julian_day_number(ymd(1970,01,01));
			static constexpr const int64_t UnixSecond = 86400LL * UnixDay;///<Calc::julian_day_number(ymd(1970,01,01));

			// Snapped from boost gregorian_calendar.ipp
			static inline unsigned long day_number(const YMDhms& ymd) noexcept {
				unsigned short a = static_cast<unsigned short>((14 - ymd.month) / 12);
				unsigned short y = static_cast<unsigned short>(ymd.year + 4800 - a);
				unsigned short m = static_cast<unsigned short>(ymd.month + 12 * a - 3);
				unsigned long  d = ymd.day + ((153 * m + 2) / 5) + 365 * y + (y / 4) - (y / 100) + (y / 400) - 32045;
				return d;
			}

			static inline YMDhms from_day_number(unsigned long dayNumber) noexcept {
				int a = dayNumber + 32044;
				int b = (4 * a + 3) / 146097;
				int c = a - ((146097 * b) / 4);
				int d = (4 * c + 3) / 1461;
				int e = c - (1461 * d) / 4;
				int m = (5 * e + 2) / 153;
				unsigned short day = static_cast<unsigned short>(e - ((153 * m + 2) / 5) + 1);
				unsigned short month = static_cast<unsigned short>(m + 3 - 12 * (m / 10));
				int year = static_cast<unsigned short>(100 * b + d - 4800 + (m / 10));
				return YMDhms(year, month, day);
			}

            static inline int64_t day_number(utctime t) noexcept {
                return (int64_t)((UnixSecond + std::chrono::duration_cast<std::chrono::seconds>(t).count()) / std::chrono::duration_cast<std::chrono::seconds>(DAY).count());
            }
            static inline utctimespan hms_seconds(int h, int m, int s,int us=0) noexcept { return deltahours(h) + deltaminutes(m) + seconds(s)+utctime(us); }



			time_zone::tz_info_t_ tz_info;
			/**\brief returns tz_info (helper for boost python really) */
			time_zone::tz_info_t_ get_tz_info() const {return tz_info;}
			/**\brief construct a timezone with standard offset, no dst, name= UTC+01 etc. */
			explicit calendar(utctimespan tz): tz_info(new time_zone::tz_info_t(tz)) {}
			explicit calendar(int tz_s=0) :tz_info(new time_zone::tz_info_t(utctimespan{seconds(tz_s)})) {}
			/**\brief construct a timezone from tz_info shared ptr provided from typically time_zone db */
			explicit calendar(time_zone::tz_info_t_ tz_info):tz_info(tz_info) {}
            calendar(calendar const&o) :tz_info(o.tz_info) {}
            calendar(calendar&&o) :tz_info(std::move(o.tz_info)) {}
            calendar& operator=(calendar const &o) {
                if (this != &o)
                    tz_info = o.tz_info;
                return *this;
            }
            calendar& operator=(calendar&&o) {
                tz_info = std::move(o.tz_info);
                return *this;
            }
            /**\brief construct a timezone based on region id
             * uses internal tz_info_database to lookup the name.
             * \param region_id like Europe/Oslo, \sa time_zone::tz_info_database
             */
            explicit calendar(std::string region_id);
            /**\brief get list of available time zone region */
            static std::vector<std::string> region_id_list();

			/**\brief construct utctime from calendar coordinates YMDhms
			 *
			 * If the YMDhms is invalid, runtime_error is thrown.
			 * Currently just trivial checks is done.
			 *
			 * \param c YMDhms that has to be valid calendar coordinates.
			 * \note special values of YMDhms, like max,min,null is mapped to corresponding utctime concepts.
			 * \sa YMDhms
			 * \return utctime
			 */
			utctime time(YMDhms c) const;
            utctime time(YWdhms c) const;

            ///<short hand for calendar::time(YMDhms)
            utctime time(int Y,int M=1,int D=1,int h=0,int m=0,int s=0,int us=0) const {
                return time(YMDhms(Y,M,D,h,m,s,us));
            }
            utctime time_from_week(int Y, int W = 1, int wd = 1, int h = 0, int m = 0, int s = 0,int us=0) const;
            /**\brief returns *utc_year* of t \note for internal dst calculations only */
			static inline int utc_year(utctime t) {
                if(t == no_utctime  ) throw std::runtime_error("year of no_utctime");
                if(t == max_utctime ) return YMDhms::YEAR_MAX;
                if(t == min_utctime ) return YMDhms::YEAR_MIN;
				return from_day_number(day_number(t)).year;
			}
            /**\brief return the calendar units of t taking timezone and dst into account
             *
             * Special utctime values, no_utctime max_utctime, min_utctime are mapped to corresponding
             * YMDhms special values.
             * \sa YMDhms
             * \return calendar units YMDhms
             */
			YMDhms  calendar_units(utctime t) const ;

            /**\brief return the calendar iso week units of t taking timezone and dst into account
            *
            * Special utctime values, no_utctime max_utctime, min_utctime are mapped to corresponding
            * YWdhms special values.
            * \sa YWdhms
            * \return calendar iso week units YWdhms
            */
            YWdhms calendar_week_units(utctime t) const;

			///< returns  0=sunday, 1= monday ... 6=sat.
			int day_of_week(utctime t) const ;

			///< returns day_of year, 1..365..
			size_t day_of_year(utctime t) const ;

			///< returns the month of t, 1..12, -1 of not valid time
			int month(utctime t) const ;

            ///< returns quarter of t, 1..4, -1 no valid time
            int quarter(utctime t) const;
            ///< returns a readable iso standard string
			std::string to_string(utctime t) const;

			///< returns a readable period with time as for calendar::to_string
			std::string to_string(utcperiod p) const;

			// calendar arithmetic
			/**\brief round down (floor) to nearest utctime with deltaT resolution.
			 *
             * If delta T is calendar::DAY,WEEK,MONTH,YEAR it do a time-zone semantically
             * correct rounding.
             * if delta T is any other number, e.g. minute/hour, the result is similar to
             * integer truncation at the level of delta T.
             * \param t utctime to be trimmed
             * \param deltaT utctimespan specifying the resolution, use calendar::DAY,WEEK,MONTH,YEAR specify calendar specific resolutions
             * \return a trimmed utctime
             */
			utctime trim(utctime t, utctimespan deltaT) const ;

			/**\brief calendar semantic add
			 *
			 *  conceptually this is similar to t + deltaT*n
			 *  but with deltaT equal to calendar::DAY,WEEK,MONTH,YEAR
			 *  and/or with dst enabled time-zone the variation of length due to dst
			 *  or month/year length is taken into account
			 *  e.g. add one day, and calendar have dst, could give 23,24 or 25 hours due to dst.
			 *  similar for week or any other time steps.
			 *
			 *  \sa calendar::diff_units
			 *
			 *  \note DST -if the calendar include dst, following rules applies:
			 *   -# transition hour 1st hour after dst has changed
			 *   -# if t and resulting t have different utc-offsets, the result is adjusted dst adjustment with the difference.
			 *
			 *
			 * \param t utctime to add n deltaT from
			 * \param deltaT that can be any, but with calendar::DAY,WEEK,MONTH,YEAR calendar semantics applies
			 * \param n number of delta T to add, can be negative
			 * \return new calculated utctime
			 */
			utctime add(utctime t, utctimespan deltaT, int64_t n) const ;

			/**\brief calculate the distance t1..t2 in specified units
			 *
			 * The function takes calendar semantics when deltaT is calendar::DAY,WEEK,MONTH,YEAR,
			 * and in addition also dst.
			 * e.g. the diff_units of calendar::DAY over summer->winter shift is 1, remainder is 0,
			 * even if the number of hours during those days are 23 and 25 summer and winter transition respectively
			 *
			 * \sa calendar::add
			 *
			 * \return (t2-t1)/deltaT, and remainder, where deltaT could be calendar units DAY,WEEK,MONTH,YEAR
			 */
			int64_t diff_units(utctime t1, utctime t2, utctimespan deltaT, utctimespan &remainder) const ;
			///< diff_units discarding remainder, \sa diff_units
			int64_t diff_units(utctime t1, utctime t2, utctimespan deltaT) const {
			    utctimespan ignore;
			    return diff_units(t1,t2,deltaT,ignore);
			}
			x_serialize_decl();
		};

        namespace time_zone {
            inline utctimespan tz_table::dst_offset(utctime t) const {
                if(!is_dst()) return  utctimespan(0);
                auto year=calendar::utc_year(t);
                if(year-start_year >= (int) dst.size()) return utctimespan(0);
                auto s=dst_start(year);
                auto e=dst_end(year);
                return (s<e? (t>=s&&t<e):(t<e || t>=s))?dt[year-start_year]:utctimespan(0);
            }
        }

	}
}
//-- serialization support: expose class keys
namespace boost {
namespace archive {
namespace sc = shyft::core;

template<class Archive>
void load(Archive& ar, sc::utctime& tp, unsigned);

template<class Archive>
void save(Archive& ar, sc::utctime const& tp, unsigned);

template<class Archive>
void serialize(Archive & ar, sc::utctime& tp, unsigned version);
#if 0
template<class Archive>
void load(Archive& ar, sc::utctimespan& tp, unsigned) ;

template<class Archive>
void save(Archive& ar, sc::utctimespan const& tp, unsigned);

template<class Archive>
void serialize(Archive & ar, sc::utctimespan& tp, unsigned version);
#endif
}
}


//x_serialize_export_key_nt(shyft::core::utctime);
//x_serialize_export_key_nt(shyft::core::utctimespan);
x_serialize_export_key_nt(shyft::core::utcperiod);
x_serialize_export_key(shyft::core::time_zone::tz_info_t);
x_serialize_export_key(shyft::core::time_zone::tz_table);
x_serialize_export_key(shyft::core::calendar);

