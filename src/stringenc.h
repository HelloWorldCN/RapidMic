/*
 * stringenc.cpp
 * 2012-04-24
 * some useful function for string manipulations,
 * including:
 *         split,            // split string by delim(such as " .,/")
 *         int2str,        // convert int to string
 *         float2str,        // convert double to string
 *         str2int,        // convert string to int
 *         str2float,        // convert string to double
 *         strtoupper,        // all to upper case
 *         strtolower,        // all to lower case
 *         //strtocapital,        // capital first character
 *
 */
#include <algorithm>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

/**
 * @brief split a string by delim
 *
 * @param str string to be splited
 * @param c delimiter, const char*, just like " .,/", white space, dot, comma,
 * splash
 *
 * @return a string vector saved all the splited world
 */
vector<string> split(string &str, const char *c);

/**
 * @brief convert a integer into string through stringstream
 *
 * @param n a integer
 *
 * @return the string form of n
 */
string int2str(int n);

string float2str(double f);

/**
 * @brief convert something to string form through stringstream
 *
 * @tparam Type Type can be int,float,double
 * @param a
 *
 * @return the string form of param a
 */
template <class Type> string tostring(Type a);

/**
 * @brief convert string to int by atoi
 *
 * @param s string
 *
 * @return the integer result
 */
int str2int(string &s);
double str2float(string &s);

/**
 * @brief do string convert through stringstream from FromType to ToType
 *
 * @tparam ToType target type
 * @tparam FromType source type
 * @param t to be converted param
 *
 * @return the target form of param t
 */
template <class ToType, class FromType> ToType strconvert(FromType t);

/**
 * @brief convert string to upper case throught transform method, also can use
 * transform method directly
 *
 * @param s
 *
 * @return the upper case result saved still in s
 */
string &strtoupper(string &s);

/**
 * @brief convert string to upper case through toupper, which transform a char
 * into upper case
 *
 * @param s
 *
 * @return the upper case result string
 */
string strtoupper(string s);

/**
 * @brief convert string to lower case throught transform method, also can use
 * transform method directly
 *
 * @param s
 *
 * @return the lower case result saved still in s
 */
string &strtolower(string &s);

/**
 * @brief convert string to lower case through tolower, which transform a char
 * into lower case
 *
 * @param s
 *
 * @return the lower case result string
 */
string strtolower(string s);
