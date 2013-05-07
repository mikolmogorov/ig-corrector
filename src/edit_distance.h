#include <vector>

template<class T>
unsigned int edit_distance(const T &s1, const T & s2) 
{
	const size_t len1 = s1.size(); 
	const size_t len2 = s2.size();
	std::vector<unsigned int> col(len2 + 1);
	std::vector<unsigned int> prevCol(len2 + 1);

	for (unsigned int i = 0; i < prevCol.size(); ++i)
	{
		prevCol[i] = i;
	}
	for (unsigned int i = 0; i < len1; ++i) 
	{
		col[0] = i + 1;
		for (unsigned int j = 0; j < len2; ++j)
		{
			unsigned int delta = s1[i]==s2[j] ? 0 : 1;
			col[j+1] = std::min( std::min(1 + col[j], 1 + prevCol[1 + j]), prevCol[j] + delta );
		}
		col.swap(prevCol);
	}
	return prevCol[len2];
}
