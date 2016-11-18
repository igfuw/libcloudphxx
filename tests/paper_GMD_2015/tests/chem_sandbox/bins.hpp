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
  // wet radius bins: .001 ... .01 ... .1 mm (25 bins in total)
  for (int i = 0; i < 25; ++i)
    ret.push_back(1e-6 * pow(10, -3 + i * .2) * si::metres); 
  return ret;
}

// focus plot locations
int ox = 0, oy=7;
std::pair<
  std::set<std::pair<int,int>>,
  std::set<std::pair<int,int>>
> focus = {
  {   // left column
    {16+ox, 20+oy}, 
    {16+ox, 37+oy}, 
    {16+ox, 39+oy}, 
    {16+ox, 44+oy},
    {16+ox, 65+oy} 
  },{ // right column
    {58+ox, 20+oy},
    {58+ox, 37+oy}, 
    {58+ox, 39+oy}, 
    {58+ox, 44+oy},
    {58+ox, 65+oy}
  }
};
