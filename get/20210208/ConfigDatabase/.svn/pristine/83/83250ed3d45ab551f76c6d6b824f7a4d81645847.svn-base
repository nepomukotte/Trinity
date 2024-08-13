#ifndef UTILS_STREAM_REQUIRE_H_INCLUDED
#define UTILS_STREAM_REQUIRE_H_INCLUDED

/**
 * Manipulator for extraction streams: requires a character to be the next one
 * in the stream and discard it.
 * It is useful to enforce formatting conventions.
 */
struct require
{
	const char ch;
	require(char ch) : ch(ch) {}
};

/**
 * Implements the "require" manipulator: ensures the character is the next in
 * the stream and discard it. Otherwise, it fails.
 */
std::istream& operator>>(std::istream& is, const require& sk)
{
	if(is.peek() == sk.ch)
	{
		is.get();
	}
	else
	{
		is.setstate(std::ios::failbit);
	}
	return is;
}

#endif // UTILS_STREAM_REQUIRE_H_INCLUDED
