// ======================= USER SETTINGS =======================
settings.render   = 5;
settings.outformat = "png"; // "pdf", "png", "svg", "eps"

string user = substr(settings.user, 0, length(settings.user) - 1);

// SET c value from user input or from here
if(user == ""){
  user = "0";
}

real cc = (real) user;

import three;

// --- Output size ---
size(230);

// --- Input files / folders ---
string folder = "results/c_" + format("%.1e", cc) + "/";   // folder containing surface + mesh files CHANGE HERE FOR DIFFERENT RESULTS
string surfaceFileName   = folder + "refined_surface_points.csv";
string edgesFileName     = folder + "edges.txt";
string boundaryFileName  = folder + "boundary_edges.txt";
string edgeRefineName    = folder + "edge_refinement.txt";

string trainLocsFileName     = "data/train_locs_4.csv"; // training locations (seed 4)
string trainResponseFileName = "data/train_response_4.csv"; // training responses (seed 4)
string coastlinesFileName    = "geo_data/coastlines.txt";
string windFileName          = "data/wind_sub.csv";

// --- Patch resolution ---
int num_points_per_curve = 10; // number of points per curve in the surface patch do not change

// --- Toggles ---
bool drawEdges = false;
bool drawAxes  = false;
bool wind      = true; // draw wind arrows

// --- Camera / lighting ---
triple cameraPos = (3, 3, 1);
triple cameraUp  = (0, 0, 1);

pen[]   lightDiffuse = new pen[]   {gray(0.8), gray(0.75)};
triple[] lightPos    = new triple[] {(2, 2, 3), (-2, -1, 2)};

// --- Visual scales ---
triple origin = 0.8*(1, -1, -1); // bottom-left corner of the merged surface
real axisLength        = 0.2;
real pointSphereRadius = 0.006;
real windArrowScale    = 0.08;

// --- Scalar normalization ---
// Previous behavior effectively assumed s in [0,1] (min=0, max=1).
// Keep that as the default to preserve the same output.
bool scalarAlreadyNormalized = true;
real scalarMin = 0;
real scalarMax = 1;

// --- Colormap / materials ---

pen interiorEdgePen = gray + 0.2bp;
pen boundaryEdgePen = gray + 0.5bp;
pen coastlinePen    = black + linewidth(0.6bp);
pen windPen         = gray(0.3);

pen surfaceEmissive = gray(0.1);
pen surfaceSpecular = black;

// --- Europe bounding box (used for coastlines / optional wind filtering) ---
real lonMin = -40, lonMax = 50;
real latMin =  30, latMax = 70;

// ======================= END SETTINGS ========================

currentprojection = perspective(cameraPos, up=cameraUp);
currentlight = light(diffuse=lightDiffuse, position=lightPos);


// =========================== HELPERS =========================
real clamp(real x, real xmin, real xmax) {
  return max(xmin, min(x, xmax));
}

// Blue → Gray → Red
pen colormap2(real t) {
  t = clamp(t, 0, 1);
  real r, g, b;

  if (t < 0.5) {
    // Blue → Gray
    real k = 2 * t;
    r = 0.0 * (1 - k) + 0.8 * k;
    g = 0.4 * (1 - k) + 0.8 * k;
    b = 1.0 * (1 - k) + 0.8 * k;
  } else {
    // Gray → Red
    real k = 2 * (t - 0.5);
    r = 0.8 * (1 - k) + 1.0 * k;
    g = 0.8 * (1 - k) + 0.0 * k;
    b = 0.8 * (1 - k) + 0.0 * k;
  }

  return rgb(r, g, b);
}

triple[] loadTriples(string filename) {
  triple[] result;
  file f = input(filename);
  while (!eof(f)) {
    string line = f;
    string[] p = split(line);
    if (p.length >= 3)
      result.push(((real)p[0], (real)p[1], (real)p[2]));
  }
  return result;
}

int[][] loadEdgeList(string filename) {
  int[][] edges;
  file f = input(filename);
  while (!eof(f)) {
    string line = f;
    string[] p = split(line);
    if (p.length >= 2)
      edges.push(new int[] {((int)p[0]), ((int)p[1])});
  }
  return edges;
}

int[] loadFlags(string filename) {
  int[] flags;
  file f = input(filename);
  while (!eof(f)) {
    string line = f;
    string[] p = split(line);
    if (p.length >= 1)
      flags.push((int)p[0]);
  }
  return flags;
}

triple geoToSphere(real lon, real lat) {
  real phi = radians(lon);
  real theta = radians(90 - lat); // colatitude
  return (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
}

// Converts from 3D point on unit sphere to (lon, lat) in degrees
pair toLonLat(triple p) {
  real lon = degrees(atan2(p.y, p.x));
  real lat = degrees(asin(p.z / length(p)));
  return (lon, lat);
}

material surfMat(pen diffusepen) {
  return material(
    diffusepen  = diffusepen,
    emissivepen = surfaceEmissive,
    specularpen = surfaceSpecular
  );
}

// =========================== AXES ============================
if (drawAxes) {
  draw(origin -- (origin + (axisLength,0,0)), Arrow3(5bp));
  label("$x$", origin + (axisLength+0.04,0,0), fontsize(9pt));

  draw(origin -- (origin + (0,axisLength,0)), Arrow3(5bp));
  label("$y$", origin + (0,axisLength+0.04,0), fontsize(9pt));

  draw(origin -- (origin + (0,0,axisLength)), Arrow3(5bp));
  label("$z$", origin + (0,0,axisLength+0.04), fontsize(9pt));
}


// ===================== LOAD EDGE METADATA ====================
int[]   bflags      = loadFlags(boundaryFileName);
triple[] nurbs_edges = loadTriples(edgeRefineName);
int[][] edges       = loadEdgeList(edgesFileName);

int num_edges = edges.length;

if (drawEdges) {
  int num_points_per_edge = num_points_per_curve;
  for (int e = 0; e < num_edges; ++e) {
    triple[] curve;
    int start = e * num_points_per_edge;
    for (int j = 0; j < num_points_per_edge; ++j)
      curve.push(nurbs_edges[start + j]);

    pen edgePen = (bflags.length > e && bflags[e] == 1) ? boundaryEdgePen : interiorEdgePen;

    for (int j = 0; j < curve.length - 1; ++j)
      draw(curve[j] -- curve[j + 1], edgePen);
  }
}


// =================== LOAD & PLOT SURFACE =====================
file surfFile = input(surfaceFileName);

int N = num_points_per_curve - 1;

triple[][][] grid;       // [cell][i][j]
real[][][]   scalarGrid; // scalar values at each point

int last_cid  = -1;
int cid_index = -1;

// Optional auto-range (disabled by default to preserve old output)
real minScalar =  1e9;
real maxScalar = -1e9;

while (!eof(surfFile)) {
  string line = surfFile;
  string[] p = split(line, ",");

  if (p.length >= 7) {
    int cid  = (int)p[0];
    int i    = (int)p[1];
    int j    = (int)p[2];
    real x   = (real)p[3];
    real y   = (real)p[4];
    real z   = (real)p[5];
    real s   = (real)p[6];

    if (cid != last_cid) {
      grid.push(new triple[N+1][N+1]);
      scalarGrid.push(new real[N+1][N+1]);
      cid_index += 1;
      last_cid = cid;
    }

    grid[cid_index][i][j] = (x, y, z);
    scalarGrid[cid_index][i][j] = s;

    if (!scalarAlreadyNormalized) {
      if (s < minScalar) minScalar = s;
      if (s > maxScalar) maxScalar = s;
    }
  }
}

if (scalarAlreadyNormalized) {
  minScalar = scalarMin;
  maxScalar = scalarMax;
}

write("Min value (blue): " + string(minScalar));
write("Max value (red): "  + string(maxScalar));


for (int c = 0; c < grid.length; ++c) {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      triple p1 = grid[c][i][j];
      triple p2 = grid[c][i+1][j];
      triple p3 = grid[c][i+1][j+1];
      triple p4 = grid[c][i][j+1];

      real s1 = scalarGrid[c][i][j];
      real s2 = scalarGrid[c][i+1][j];
      real s3 = scalarGrid[c][i+1][j+1];
      real s4 = scalarGrid[c][i][j+1];

      // Triangle 1: p1-p2-p3
      real t1 = ((s1 + s2 + s3) / 3 - minScalar) / (maxScalar - minScalar + 1e-10);
      pen color1 = colormap2(round(t1*1000)/1000); // round to 3 decimal places
      draw(surface(p1--p2--p3--cycle), surfacepen = surfMat(color1));

      // Triangle 2: p1-p3-p4
      real t2 = ((s1 + s3 + s4) / 3 - minScalar) / (maxScalar - minScalar + 1e-10);
      pen color2 = colormap2(round(t2*1000)/1000); // round to 3 decimal places
      draw(surface(p1--p3--p4--cycle), surfacepen = surfMat(color2));
    }
  }
}


// ===================== TRAINING POINTS =======================
file locFile      = input(trainLocsFileName);
file responseFile = input(trainResponseFileName);

bool skipHeader = true;
while (!eof(locFile) && !eof(responseFile)) {
  string locLine  = locFile;
  string respLine = responseFile;

  if (skipHeader) { skipHeader = false; continue; }

  string[] locParts  = split(locLine, ",");
  string[] respParts = split(respLine, ",");
  if (locParts.length < 4 || respParts.length < 2) continue;

  real x = (real)locParts[1];
  real y = (real)locParts[2];
  real z = (real)locParts[3];
  real s = (real)respParts[1];

  real t = (s - minScalar) / (maxScalar - minScalar + 1e-10);
  pen color = colormap2(t);

  triple pt = (x, y, z);
  draw(shift(pt)*scale3(pointSphereRadius)*unitsphere,
       surfacepen = material(diffusepen=color, emissivepen=black, specularpen=black));
}


// ======================= COASTLINES ==========================
file coastFile = input(coastlinesFileName);

triple lastpt;
bool first = true;

while (!eof(coastFile)) {
  string line = coastFile;
  line = replace(line, "\r", ""); // handle CRLF

  if (line == "") {
    first = true; // new polygon
    continue;
  }

  string[] parts = split(line);
  if (parts.length < 2) continue;

  real lon = (real)parts[0];
  real lat = (real)parts[1];

  // Keep only points within Europe box
  if (lon < lonMin || lon > lonMax || lat < latMin || lat > latMax) {
    first = true;
    continue;
  }

  triple pt = geoToSphere(lon, lat);
  if (!first) draw(lastpt -- pt, coastlinePen);
  lastpt = pt;
  first = false;
}


// ========================= WIND ==============================
file windFile = input(windFileName);
bool skipWindHeader = true;

while (!eof(windFile) && wind) {
  string line = windFile;
  if (skipWindHeader) { skipWindHeader = false; continue; }

  string[] parts = split(line, ",");
  if (parts.length < 6) continue;

  triple p = ((real)parts[0], (real)parts[1], (real)parts[2]);
  triple b = ((real)parts[3], (real)parts[4], (real)parts[5]);

  // Optional geographic filter (disabled previously; keep it optional here)
  //pair lonlat = toLonLat(p);
  //real lon = lonlat.x;
  //real lat = lonlat.y;
  //if (lon < lonMin || lon > lonMax || lat < latMin || lat > latMax) continue;

  draw(p -- (p + windArrowScale * b), windPen, Arrow3(3bp));
}
