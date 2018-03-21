#include "frame.h"
#include "system.h"
#include "atom.h"
#include "radialDistribution.h"
#include <vector>
#include <iostream>

using namespace std;

static FILE *fp = NULL;

static void writeString(const char *str)
{
    fprintf(fp, "%s",str);
}

static void writeFloat(float val)
{
  char str[128];
  sprintf(str, "%20.12e ", val);
  fprintf(fp, "%s",str);
}

RDF::RDF()
{
};

RDF::RDF(System& a_system)
{
  m_system = a_system;
  m_rdf.resize(m_system.getNumPairs());
  m_pairCounter.resize(m_system.getNumPairs());
  m_maxDist = 14.0;
  m_numBins = 500;
  for (int i=0; i<m_system.getNumPairs(); i++)
    {
      m_rdf[i].resize(m_numBins);
    }
  m_binSize = m_maxDist / m_numBins;
};

void RDF::sample(Frame& a_frame)
{
  double pairDistance;
  for (int pairIdx=0; pairIdx<m_system.getNumPairs(); pairIdx++) {
    pair<unsigned int, unsigned int > pair = m_system.getPairCorrelation(pairIdx);
    for (int atomTypeI=0; atomTypeI<m_system.getNumOfType(pair.first); atomTypeI++)
      {
	for (int atomTypeJ=0; atomTypeJ<m_system.getNumOfType(pair.second); atomTypeJ++)
	  {
	    pairDistance = a_frame.computeDistance(m_system.getIndexOfType(pair.first,atomTypeI), m_system.getIndexOfType(pair.second,atomTypeJ));
	    if (pairDistance < m_maxDist)
	      {
		binPairDistance(pairIdx, pairDistance);
		incrementCounter(pairIdx);
	      }
	  }
      }
  }
}
void RDF::binPairDistance(unsigned int a_pair, double a_distance)
{
  int bin = floor(a_distance / m_binSize);
  m_rdf[a_pair][bin]++;
}

void RDF::incrementCounter(unsigned int a_pair)
{
  m_pairCounter[a_pair]++;
}

void RDF::normalize()
{
  for (int i=0; i<m_system.getNumPairs(); i++)
    {
      for (int j=0; j<m_numBins; j++)
	{
	  //m_rdf[i][j] /= m_pairCounter[i];
	  m_rdf[i][j] /= (4./3.)*M_PI*(pow(j+1,3)-pow(j,3))*pow(m_binSize,3);
	}
    }
}

void RDF::print()
{
  {
    for (int i=0; i<m_system.getNumPairs(); i++)
      {
	for (int j=0; j<m_numBins; j++)
	  {
	    cout << j*m_binSize << " " << m_rdf[i][j] << endl;
	  }
	cout << endl << endl;
      }
  }

}

const unsigned int RDF::getNumPairs() const
{
  return m_system.getNumPairs();
}
const unsigned int RDF::getNumBins() const
{
  return m_numBins;
}

const double RDF::getBinSize() const
{
  return m_binSize;
}

const double RDF::getRDFElement(int i, int j) const
{
  return m_rdf[i][j];
}

const char* RDFWrite(RDF* a_rdf, const char* a_filename)
{
  char full_filename[1024];
  if (strstr(a_filename, ".out") != NULL)
    {
      strcpy(full_filename, a_filename);
    }
  else
    {
      sprintf(full_filename, "%s.out", a_filename);
    }
  fp = fopen(full_filename, "w+");

  for (unsigned int iBin=0; iBin<a_rdf->getNumBins(); iBin++)
    {
      char str[128];
      sprintf(str, "%f", iBin*a_rdf->getBinSize());
      writeString(str);
      for (unsigned int jPair=0; jPair<a_rdf->getNumPairs(); jPair++)
	{
	  char str[128];
	  sprintf(str, " %f", a_rdf->getRDFElement(jPair,iBin));
	  writeString(str);
	}
      writeString("\n");
    }
  fclose(fp);
  return a_filename;

}

