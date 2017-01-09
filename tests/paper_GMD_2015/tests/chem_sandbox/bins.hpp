#pragma once

// bin sizes for calc and plot
vector<quantity<si::length>> bins_dry()
{
  vector<quantity<si::length>> ret;
  // dry radius bins: .001 ... .01 ... 10 (40 bins in total)
  for (int i = 0; i < 40; ++i)
    ret.push_back(1e-6 * pow(10, -3 + i * .1) * si::metres);
  return ret;
}

vector<quantity<si::length>> bins_wet()
{
  vector<quantity<si::length>> ret;
  // wet radius bins: .001um ... .01 ... .1 mm + 5 more bins (55 bins in total)
  for (int i = 0; i < 55; ++i)
    ret.push_back(1e-6 * pow(10, -3 + i * .1) * si::metres); 
  return ret;
}

// focus plot locations
int ox = 0, oy=0;
std::pair<
  std::set<std::pair<int,int>>,
  std::set<std::pair<int,int>>
> focus = {
  {   // left column
    {12+ox, 37+oy},
    {12+ox, 7+oy}
  },{ // right column
    {31+ox, 7+oy},
    {41+ox, 7+oy}
  }
};
