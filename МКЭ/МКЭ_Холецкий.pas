uses Graph3D, GraphABC, NumLibABC, System, Cholesky;


var
  R := 1;
  dr := R / 10; // Максимум 68. Оценка Pn = 3.333*x^2 (грубая), для n<30: Pn = 5.4501*x^1.8705
  
  r_dt_nt := new List<(real, integer)>;
  b_arr: array of real;
  u_arr: M_Zip;
  u_arr_: Matrix;
  
  Nsum: integer;
  
  // Переменные для Draw
  Z_Ax_min := -1.0;
  Z_Ax_max := 1.0;
  IsDraw := false;
  m_GR := DiffuseMaterial(Graph3D.ARGB(10, 125, 125, 125));
  F_All := new List<TriangleT>;
  U_All := new List<TriangleT>;
  B_All := new List<TriangleT>;
  D_All := new List<SegmentsT>;
  N_All := new List<TextT>;


type
  TriDel = class
  private
  
  public
    N: integer;
    P0: point3D;
    P1: point3D;
    P2: point3D;
    U0: point3D;
    U1: point3D;
    U2: point3D;
    D0: point3D;
    D1: point3D;
    D2: point3D;
    NegTr: array[0..2] of TriDel;
    
    ///Определяет, какому кольцу принадлежит среднее ребро (true - верхнему)
    TrUp: boolean; 
    Val: real;
    Fl: boolean := false;
    constructor Create(x0, y0, x1, y1, x2, y2: real; TrUp_: boolean; N_: integer);
    begin
      P0 := new Point3D(x0, y0, 0);
      P1 := new Point3D(x1, y1, 0);
      P2 := new Point3D(x2, y2, 0);
      U0 := P0;
      U1 := P1;
      U2 := P2;
      D0 := P0;
      D1 := P1;
      D2 := P2;
      TrUp := TrUp_;
      
      N := N_;
    end;
  end;

var
  TrArr := new List<List<TriDel>>;

var
  a0_, a1_, a2_, b0_, b1_, b2_, c0_, c1_, c2_, d0_, d1_, d2_: real;

var
  T0_: TriDel;

// -------------------------------------
//-----------Триангурирование-----------
// -------------------------------------

function d_ab(xa, ya, xb, yb: real): real;
begin
  result := Sqrt(Sqr(xa - xb) + Sqr(ya - yb));
end;

function XYtoN(P: Point3D): integer;
begin
  var Rnow := Sqrt(P.X * P.X + P.Y * P.Y);
  var i := ((R - Rnow) / dr).Round;
  var k := ((R / dr - i).Round.IsOdd) ? r_dt_nt[i][0] / 2 : 0;
  var dtnow := Math.Atan2(P.Y, P.X) - k + r_dt_nt[i][0] / 10;
  if dtnow < 0 then dtnow += 2 * Pi;
  var Nnow := (dtnow / r_dt_nt[i][0]).Round;
  Result := Nnow + r_dt_nt[i][1];
end;

function Delone(R: integer; dr: real): integer;
var
  x1, y1, x2, y2, xa, ya, xb, yb: real;
begin
  
  var nt1 := (2 * Pi * R / dr).Round;
  var nt2 := nt1;
  
  var i := -1;
  var Nsum := 0;
  var NN := 0;
  
  for var i1 := (R / dr).Round downto 1 do
  begin
    var R1 := i1 * dr;
    var dt1 := (2 * Pi / nt1);
    //if (i1 mod 10 = 0) or (2 * Pi * R1 / nt1 < dr / 4 * 3) then begin
    if (2 * Pi * R1 / nt1 < dr / 4 * 3) then begin
      nt1 := (2 * Pi * R1 / dr).Round;
      dt1 := (2 * Pi / nt1);
    end;
    var kt1 := (i1.IsOdd) ? dt1 / 2 : 0;
    
    
    var i2 := i1 - 1;
    var R2 := i2 * dr;
    var dt2 := dt1;
    //if (i2 mod 10 = 0) or (2 * Pi * R2 / nt2 < dr / 4 * 3) then begin
    if (2 * Pi * R2 / nt2 < dr / 4 * 3) then begin
      nt2 := (2 * Pi * R2 / dr).Round;
      dt2 := (2 * Pi / nt2);
    end;
    var kt2 := (i2.IsOdd) ? dt2 / 2 : 0;
    
    xa := (R1 * cos(0 + kt1));
    ya := (R1 * sin(0 + kt1));
    xb := (R2 * cos(0 + kt2));
    yb := (R2 * sin(0 + kt2));
    
    
    var t1 := dt1;
    var t2 := dt2;
    
    x1 := (R1 * cos(t1 + kt1));
    y1 := (R1 * sin(t1 + kt1));
    if R2 > dr / 2 then begin
      x2 := (R2 * cos(t2 + kt2));
      y2 := (R2 * sin(t2 + kt2));
    end
    else begin
      x2 := 0;
      y2 := 0;
      xb := 0;
      yb := 0;
      dt2 := dt1;
    end;
    
    var xa0 := xa;
    var ya0 := ya;
    var xb0 := xb;
    var yb0 := yb;
    
    
    // + в массив для расчета N по XY
    // Для X=0, Y=0 записывается после цикла!!!
    r_dt_nt.Add((dt1, Nsum));
    Nsum += nt1;
    
    i += 1;
    TrArr.Add(new List<TriDel>);
    
    repeat
      if ((d_ab(xa, ya, x2, y2) < d_ab(xb, yb, x1, y1)) and (R2 > dr / 2)) then begin
        
        // Базовый треугольник
        TrArr[i].Add(new TriDel(x2, y2, xa, ya, xb, yb, false, NN));
        NN += 1;
        xb := x2;
        yb := y2;
        t2 += dt2;
        x2 := (R2 * cos(t2 + kt2));
        y2 := (R2 * sin(t2 + kt2));
        
      end 
      else begin
        
        TrArr[i].Add(new TriDel(x1, y1, xb, yb, xa, ya, true, NN));
        NN += 1;
        xa := x1;
        ya := y1;
        t1 += dt1;
        x1 := (R1 * cos(t1 + kt1));
        y1 := (R1 * sin(t1 + kt1));       
      end;
    until ((Abs(xa0 - xa) < dr / 100) and (Abs(ya0 - ya) < dr / 100) and (Abs(xb0 - xb) < dr / 100) and (Abs(yb0 - yb) < dr / 100));
  end;
  
  // + в массив для расчета N по XY для 0,0
  r_dt_nt.Add((2 * Pi, Nsum));
  Result := Nsum + 1;
end;

// Для ввода координаты Z используется картинка в грейскале, где Цвет 90 = 0 по оси Z
// Картинка  2001х2001 пиксель!!!!! Из-за округления.
procedure ZfromPic(file_name: string);
const
  scaleZ = 1000;
  scaleXY = 1000;
begin
  var pic := new Picture(1, 1);
  pic.Load(file_name); // leon.jpg skull.jpg  t1.jpg
  Println(pic.Width, pic.Height);
  
  
  for var k := 0 to TrArr.Count - 1 do
  begin
    for var j := 0 to TrArr[k].Count - 1 do
    begin
      TrArr[k][j].P0.Z := (pic.GetPixel((TrArr[k][j].P0.X * scaleXY).Round + scaleXY, (TrArr[k][j].P0.Y * scaleXY).Round + scaleXY).R - 90) / scaleZ;
      TrArr[k][j].P1.Z := (pic.GetPixel((TrArr[k][j].P1.X * scaleXY).Round + scaleXY, (TrArr[k][j].P1.Y * scaleXY).Round + scaleXY).R - 90) / scaleZ;
      TrArr[k][j].P2.Z := (pic.GetPixel((TrArr[k][j].P2.X * scaleXY).Round + scaleXY, (TrArr[k][j].P2.Y * scaleXY).Round + scaleXY).R - 90) / scaleZ;
      
      var Zmin, Zmax: real;
      Zmin := Min(TrArr[k][j].P0.Z, TrArr[k][j].P1.Z, TrArr[k][j].P2.Z);
      Zmax := Max(TrArr[k][j].P0.Z, TrArr[k][j].P1.Z, TrArr[k][j].P2.Z);
      if Z_Ax_min > Zmin then Z_Ax_min := Zmin;
      if Z_Ax_max < Zmax then Z_Ax_max := Zmax;  
    end;
  end;
end;

procedure ZfromFunc;
  function user_f(P0: Point3D): real;
  begin
    var X := P0.X / R;
    var Y := P0.Y / R;
    //Result := ((X+0) * (X+0) + Y * Y) - 5;
    //Result := ((X + 0.5) * (X + 0.5) + Y * Y) - 5;
    //Result := ((X + 0.5) * (X + 0.5) + Y * Y) - 5;
    //Result := Sin(Pi * (X * X + Y * Y)) + 4 * Pi * Pi * X * X * Sin(Pi * (X * X + Y * Y)) - 4 * Pi * Cos(Pi * (X * X + Y * Y)) + 4 * Pi * Pi * Y * Y * Sin(Pi * (X * X + Y * Y));  
    //Result := 1+Sin(6*(X*X+Y*Y));
    Result := 100*(X*X+Y*Y) * (1 + cos(5*arctan(Y/(X+0.000001))));
    //Result := (X+Y)/2;
    //Result := (X * X + Y * Y) - 8*X * X * X * X;
    //Result := 2*Sin(X/Pi/2);
    //Result := -5;
    //Result := -1+((X + 0.0) * (X + 0.0) + Y * Y);
    //Result := Sqrt(1-((X + 0.0) * (X + 0.0) + Y * Y));
    
    //=========== для -u''(X,Y), u = sin(pi(x*x+y*y)) ==================
    //Result := -4*Pi*Pi*X*X*Sin(Pi*(X*X+Y*Y))+4*Pi*Cos(Pi*(X*X+Y*Y))-4*Pi*Pi*Y*Y*Sin(Pi*(X*X+Y*Y)); 
    
  end;

begin
  for var k := 0 to TrArr.Count - 1 do
  begin
    for var j := 0 to TrArr[k].Count - 1 do
    begin
      TrArr[k][j].P0.Z := user_f(TrArr[k][j].P0);
      TrArr[k][j].P1.Z := user_f(TrArr[k][j].P1);
      TrArr[k][j].P2.Z := user_f(TrArr[k][j].P2);
      
      var Zmin, Zmax: real;
      Zmin := Min(TrArr[k][j].P0.Z, TrArr[k][j].P1.Z, TrArr[k][j].P2.Z);
      Zmax := Max(TrArr[k][j].P0.Z, TrArr[k][j].P1.Z, TrArr[k][j].P2.Z);
      if Z_Ax_min > Zmin then Z_Ax_min := Zmin;
      if Z_Ax_max < Zmax then Z_Ax_max := Zmax;  
    end;
  end;
end;




// ===================== Рисование =============================
// =============================================================

procedure Draw3D_Ax();
begin
  // Оси XYZ - рисование
  Z_Ax_min := Z_Ax_min;
  Z_Ax_max := Z_Ax_max;
  
  var dn_Ax := 1 / 100;
  // Ось Z
  for var j := ((Z_Ax_min - 1 / 2) / dn_Ax).Round to ((Z_Ax_max + 1 / 2) / dn_Ax).Round do
    Polyline3D(Seq(new Point3D(0, 0, j * dn_Ax), new Point3D(0, 0, (j + 1) * dn_Ax)), 1, Graph3D.ARGB(200, 0, 0, 255));
  // Ось X
  for var j := ((-1 - 1 / 5) / dn_Ax).Round to ((1 + 1 / 5) / dn_Ax).Round do
    Polyline3D(Seq(new Point3D(j * dn_Ax, 0, 0), new Point3D((j + 1) * dn_Ax, 0, 0)), 1, Graph3D.ARGB(200, 255, 0, 0));
  // Ось Y
  for var j := ((-1 - 1 / 5) / dn_Ax).Round to ((1 + 1 / 5) / dn_Ax).Round do
    Polyline3D(Seq(new Point3D(0, j * dn_Ax, 0), new Point3D(0, (j + 1) * dn_Ax, 0)), 1, Graph3D.ARGB(200, 0, 255, 0));
  
  //GraphAbc.GraphABCControl.Cursor := system.windows.forms.cursors.Cross;
end;

procedure Draw3D_B();
var
  K0, K1, K2: Point3D;
  m_0 := DiffuseMaterial(Graph3D.ARGB(200, 255, 255, 0));
  m_1 := DiffuseMaterial(Graph3D.ARGB(200, 255, 0, 255));
  m_2 := DiffuseMaterial(Graph3D.ARGB(200, 0, 255, 255));
  a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2: real;

begin
  IsDraw := true;
  
  if B_All.Count = 0 then begin
    var KK0, KK1, KK2: Point3D;
    
    for var k := 0 to TrArr.Count - 1 do
      for var j := 0 to TrArr[k].Count - 1 do
      begin
        K0 := TrArr[k][j].P0;
        K1 := TrArr[k][j].P1;
        K2 := TrArr[k][j].P2;
        
        // Коэффиценты при базисе (K0-K1) вида (a0x+b0y+c0)/d0
        d0 := (K2.X - K1.X) * (K1.Y - K2.Y) - (K0.X - K1.X) * (K1.Y - K2.Y) - (K0.Y - K2.Y) * (K2.X - K1.X); 
        a0 := K2.Y - K1.Y;
        b0 := K1.X - K2.X;
        c0 := (K2.X - K1.X) * (K1.Y - K2.Y) + K1.X * (K1.Y - K2.Y) + K2.Y * (K2.X - K1.X);
        
        // Коэффиценты при базисе (K1-K2) вида (a1x+b1y+c1)/d1
        d1 := (K0.X - K2.X) * (K2.Y - K0.Y) - (K1.X - K2.X) * (K2.Y - K0.Y) - (K1.Y - K0.Y) * (K0.X - K2.X); 
        a1 := K0.Y - K2.Y;
        b1 := K2.X - K0.X;
        c1 := (K0.X - K2.X) * (K2.Y - K0.Y) + K2.X * (K2.Y - K0.Y) + K0.Y * (K0.X - K2.X);
        
        // Коэффиценты при базисе (K2-K0) вида (a2x+b2y+c2)/d2
        d2 := (K1.X - K0.X) * (K0.Y - K1.Y) - (K2.X - K0.X) * (K0.Y - K1.Y) - (K2.Y - K1.Y) * (K1.X - K0.X); 
        a2 := K1.Y - K0.Y;
        b2 := K0.X - K1.X;
        c2 := (K1.X - K0.X) * (K0.Y - K1.Y) + K0.X * (K0.Y - K1.Y) + K1.Y * (K1.X - K0.X);
        
        
        K0.Z := (a0 * K0.X + b0 * K0.Y + c0) / d0;
        K1.Z := (a0 * K1.X + b0 * K1.Y + c0) / d0;
        K2.Z := (a0 * K2.X + b0 * K2.Y + c0) / d0;
        KK0 := K0;
        KK1 := K1;
        KK2 := K2;
        B_All.Add(Triangle(KK0, KK1, KK2, m_0)); 
        
        K0.Z := (a1 * K0.X + b1 * K0.Y + c1) / d1;
        K1.Z := (a1 * K1.X + b1 * K1.Y + c1) / d1;
        K2.Z := (a1 * K2.X + b1 * K2.Y + c1) / d1;
        KK0 := K0;
        KK1 := K1;
        KK2 := K2;
        B_All.Add(Triangle(KK0, KK1, KK2, m_1)); 
        
        K0.Z := (a2 * K0.X + b2 * K0.Y + c2) / d2;
        K1.Z := (a2 * K1.X + b2 * K1.Y + c2) / d2;
        K2.Z := (a2 * K2.X + b2 * K2.Y + c2) / d2;
        KK0 := K0;
        KK1 := K1;
        KK2 := K2;
        B_All.Add(Triangle(KK0, KK1, KK2, m_2));
      end;
    
  end else
    foreach var f in B_All do f.Visible := not (f.Visible);
  
  IsDraw := false;
end;

procedure Draw3D_F();
var
  K0, K1, K2: Point3D;
  //m_ := DiffuseMaterial(Graph3D.ARGB(255, 190, 255, 190));
  m_ := RainbowMaterial;
begin
  IsDraw := true;
  
  if F_All.Count = 0 then begin
    for var k := 0 to TrArr.Count - 1 do
      for var j := 0 to TrArr[k].Count - 1 do
      begin
        //Polygon3D(Seq(TrArr[k][j].P0, TrArr[k][j].P1, TrArr[k][j].P2), 1, GrayColor(0)); // Ребра
        K0 := TrArr[k][j].P0;
        K1 := TrArr[k][j].P1;
        K2 := TrArr[k][j].P2;
        F_All.Add(Triangle(K0, K1, K2, m_));         // Треугольники
      end;
  end else
    foreach var f in F_All do f.Visible := not (f.Visible);
  
  IsDraw := false;
end;

procedure Draw3D_U();
var
  K0, K1, K2: Point3D;
  m_ := DiffuseMaterial(Graph3D.ARGB(255, 255, 180, 255));
  //m_ := RainbowMaterial;
  
  function P_P(P1, P2: Point3D; n: String): Point3D;
  var
    K: real;
  begin
    if (n = 'YZ') or (n = 'ZY') then K := (P1.X) / (P1.X - P2.X);
    if (n = 'XZ') or (n = 'ZX') then K := (P1.Y) / (P1.Y - P2.Y);
    if (n = 'XY') or (n = 'YX') then K := (P1.Z) / (P1.Z - P2.Z);
    if (Abs(K) > 1) or (K < 0) or (real.IsNaN(K)) then Result := BadPoint
    else begin
      Result := new Point3D((P2.X - P1.X) * K + P1.X, (P2.Y - P1.Y) * K + P1.Y, (P2.Z - P1.Z) * K + P1.Z);
    end;
    //Println(V1,V2,Result, K);
  end;

begin
  IsDraw := true;
  
  if U_All.Count = 0 then begin
    for var k := 0 to TrArr.Count - 1 do
      for var j := 0 to TrArr[k].Count - 1 do
      begin
        //Polygon3D(Seq(TrArr[k][j].U0, TrArr[k][j].U1, TrArr[k][j].U2), 3, GrayColor(0)); // Ребра
        K0 := TrArr[k][j].U0;
        K1 := TrArr[k][j].U1;
        K2 := TrArr[k][j].U2;
        
        var Tr_now := Triangle(K0, K1, K2, m_);
        U_All.Add(Tr_now);         // Треугольники
        
        if (Min(Tr_now.P1.Y, Tr_now.P2.Y, Tr_now.P3.Y) <= 0) and (Max(Tr_now.P1.Y, Tr_now.P2.Y, Tr_now.P3.Y) >= 0) then begin
          var G1 := P_P(Tr_now.P1, Tr_now.P2, 'XZ');
          var G2 := P_P(Tr_now.P2, Tr_now.P3, 'XZ');
          var G3 := P_P(Tr_now.P3, Tr_now.P1, 'XZ');
          if G1 = BadPoint then Polygon3D(Seq(G2, G3), 3, Graph3D.RGB(255, 0, 255))
          else if G2 = BadPoint then Polygon3D(Seq(G1, G3), 3, Graph3D.RGB(255, 0, 255))
          else if G3 = BadPoint then Polygon3D(Seq(G1, G2), 3, Graph3D.RGB(255, 0, 255))
          else if G1 <> G2 then Polygon3D(Seq(G1, G2), 3, Graph3D.RGB(255, 0, 255))
          else Polygon3D(Seq(G2, G3), 3, Graph3D.RGB(255, 0, 255));
        end;
        
        if (Min(Tr_now.P1.X, Tr_now.P2.X, Tr_now.P3.X) <= 0) and (Max(Tr_now.P1.X, Tr_now.P2.X, Tr_now.P3.X) >= 0) then begin
          var G1 := P_P(Tr_now.P1, Tr_now.P2, 'YZ');
          var G2 := P_P(Tr_now.P2, Tr_now.P3, 'YZ');
          var G3 := P_P(Tr_now.P3, Tr_now.P1, 'YZ');
          if G1 = BadPoint then Polygon3D(Seq(G2, G3), 3, Graph3D.RGB(0, 255, 255))
          else if G2 = BadPoint then Polygon3D(Seq(G1, G3), 3, Graph3D.RGB(0, 255, 255))
          else if G3 = BadPoint then Polygon3D(Seq(G1, G2), 3, Graph3D.RGB(0, 255, 255))
          else if G1 <> G2 then Polygon3D(Seq(G1, G2), 3, Graph3D.RGB(0, 255, 255))
          else Polygon3D(Seq(G2, G3), 3, Graph3D.RGB(0, 255, 255));
        end;
        
      end;
  end else
    foreach var f in U_All do f.Visible := not (f.Visible);
  
  IsDraw := false;
end;


// Рисуем триагуляцию Делоне 
procedure Draw3D_D();
var
  K0, K1, K2: Point3D;
begin
  IsDraw := true;
  
  if D_All.Count = 0 then begin
    for var k := 0 to TrArr.Count - 1 do
      for var j := 0 to TrArr[k].Count - 1 do
      begin
        var PP0 := TrArr[k][j].P0;
        PP0.Z := 0;
        var PP1 := TrArr[k][j].P1;
        PP1.Z := 0;
        var PP2 := TrArr[k][j].P2;
        PP2.Z := 0;
        
        D_All.Add(Polygon3D(Seq(PP0, PP1, PP2), 1, Graph3D.ARGB(255, 0, 0, 255)));
        
      end;
  end else begin
    foreach var f in D_All do f.Destroy;
    D_All.Clear();
  end;
  
  IsDraw := false;
end;

// Нумерация вершин без повторов (не тест)
// 
procedure Draw3D_N();
var
  K0, K1, K2: Point3D;
begin
  IsDraw := true;
  
  var un_arr := ArrFill(Nsum * Nsum, true);
  
  if N_All.Count = 0 then begin
    for var k := 0 to TrArr.Count - 1 do
      for var j := 0 to TrArr[k].Count - 1 do
      begin
        var PP0 := TrArr[k][j].P0;
        PP0.Z := 0;
        var PP1 := TrArr[k][j].P1;
        PP1.Z := 0;
        var PP2 := TrArr[k][j].P2;
        PP2.Z := 0;
        
        N_All.Add(Text3D(new Point3D((PP0.X + PP1.X + PP2.X) / 3, (PP0.Y + PP1.Y + PP2.Y) / 3, 0), TrArr[k][j].N.ToString, (dr / 4)));
        
        if un_arr[XYtoN(PP0)] then begin
          N_All.Add(Text3D(PP0, XYtoN(PP0).ToString, (dr / 4)));
          un_arr[XYtoN(PP0)] := false;
        end;
        if un_arr[XYtoN(PP1)] then begin
          N_All.Add(Text3D(PP1, XYtoN(PP1).ToString, (dr / 4)));
          un_arr[XYtoN(PP1)] := false;
        end;
        if un_arr[XYtoN(PP2)] then begin
          N_All.Add(Text3D(PP2, XYtoN(PP2).ToString, (dr / 4)));
          un_arr[XYtoN(PP2)] := false;
        end;
        
        //        N_All.Add(BillboardText(new Point3D((PP0.X + PP1.X + PP2.X) / 3, (PP0.Y + PP1.Y + PP2.Y) / 3, 0), TrArr[k][j].N.ToString, 15));
        //        N_All.Add(BillboardText(PP0, XYtoN(PP0).ToString, 15));
        //        N_All.Add(BillboardText(PP1, XYtoN(PP1).ToString, 15));
        //        N_All.Add(BillboardText(PP2, XYtoN(PP2).ToString, 15));
      end;
  end else begin
    foreach var f in N_All do f.Destroy;
    N_All.Clear();
  end;
  
  IsDraw := false;
end;


// Тестовая нумерация вершин 
//(долго, т.к. нумеруются повторно для проверки сбоев нумерации и дубликатов - текст накладывается)
procedure Draw3D_TN();
var
  K0, K1, K2: Point3D;
begin
  IsDraw := true;
  
  if N_All.Count = 0 then begin
    for var k := 0 to TrArr.Count - 1 do
      for var j := 0 to TrArr[k].Count - 1 do
      begin
        var PP0 := TrArr[k][j].P0;
        PP0.Z := 0;
        var PP1 := TrArr[k][j].P1;
        PP1.Z := 0;
        var PP2 := TrArr[k][j].P2;
        PP2.Z := 0;
        
        N_All.Add(Text3D(new Point3D((PP0.X + PP1.X + PP2.X) / 3 , (PP0.Y + PP1.Y + PP2.Y) / 3 , 0), TrArr[k][j].N.ToString, (dr / 4 )));
        N_All.Add(Text3D(PP0 , XYtoN(PP0).ToString, (dr / 4 )));
        N_All.Add(Text3D(PP1 , XYtoN(PP1).ToString, (dr / 4 )));
        N_All.Add(Text3D(PP2 , XYtoN(PP2).ToString, (dr / 4 )));
        
        //        N_All.Add(BillboardText(new Point3D((PP0.X + PP1.X + PP2.X) / 3 , (PP0.Y + PP1.Y + PP2.Y) / 3 , 0), TrArr[k][j].N.ToString, 15));
        //        N_All.Add(BillboardText(PP0 , XYtoN(PP0).ToString, 15));
        //        N_All.Add(BillboardText(PP1 , XYtoN(PP1).ToString, 15));
        //        N_All.Add(BillboardText(PP2 , XYtoN(PP2).ToString, 15));
      end;
  end else begin
    foreach var f in N_All do f.Destroy;
    N_All.Clear();
  end;
  
  IsDraw := false;
end;



// -------------------------------------
// -----------Интегрирование------------
// -------------------------------------
function Fx(N0, N1: Point3D) := new Polynom(N0.Y - (N1.Y - N0.Y) / (N1.X - N0.X) * N0.X, (N1.Y - N0.Y) / (N1.X - N0.X)); // y=ax+b


procedure integral_triangle(T0: TriDel); // Помещает коэффиценты в матрицу при треугольнике
var
  a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2, koef_u, koef_b: real;
  K0, K1, K2, KMin, KMid, KMax, KNew: Point3D;
  FxLU, FxLD, FxRU, FxRD: Polynom;
  N0, N1, N2: integer;
  
  function integrate_basis(ai, aj, bi, bj, ci, cj, di, dj: real; diff: boolean): real;
  var
    res, temp: real;
  begin
    res := 0;
    temp := 0;
    var FxRes0 := new Polynom(ci * bj + cj * bi, ai * bj + aj * bi);
    var FxRes1 := new Polynom(ci, ai);
    var FxRes2 := new Polynom(cj, aj);
    var diff_koef := 0.0;
    if diff then diff_koef := (ai * aj + bi * bj) / di / dj; // Векторное
    
    var FxRes: Polynom;
    var d_NaN := dr / 1000;
    
    // ============ Отладка u-u'', u, -u'' ==============
    
    var K_U := 1;
    var K_dU := 1;
    
    // 1. u-u''; 2. u ; 3. -u'' 
    // Не забывать менять F(x)!!!!!
    var cc := 1;
    case cc of
      1:  ; // 1. u-u''
      2: K_dU := 0; // 2. u
      3: begin K_dU := -1; if diff then K_U := 0; end; // 3. -u'' 
    end;
    // ===================================================
    
    
    if Abs(KMid.X - KMin.X) > d_NaN then begin
      FxRes := K_U * ((bi * bj) / 3 * (FxLU * FxLU * FxLU - FxLD * FxLD * FxLD) + (FxRes0) / 2 * (FxLU * FxLU - FxLD * FxLD) + (FxRes1 * FxRes2) * (FxLU - FxLD)) / di / dj + K_dU * diff_koef * (FxLU - FxLD);
      
      FxRes := FxRes.PInt;
      //FxRes.PrintlnBeauty;
      res += FxRes.Value(KMid.X) - FxRes.Value(KMin.X);
    end;
    
    
    
    if Abs(KMax.X - KMid.X) > d_NaN then begin
      FxRes := K_U * ((bi * bj) / 3 * (FxRU * FxRU * FxRU - FxRD * FxRD * FxRD) + (FxRes0) / 2 * (FxRU * FxRU - FxRD * FxRD) + (FxRes1 * FxRes2) * (FxRU - FxRD)) / di / dj + K_dU * diff_koef * (FxRU - FxRD);
      FxRes := FxRes.PInt;
      //FxRes.PrintlnBeauty;
      res += FxRes.Value(KMax.X) - FxRes.Value(KMid.X);
    end;
    
    
    Result := res;
  end;

begin
  K0 := T0.P0;
  K1 := T0.P1;
  K2 := T0.P2;
  
  KMin := K0;
  KMid := K1;
  KMax := K2;
  if (KMin.X > KMid.X) then Swap(KMin, KMid);
  if (KMin.X > KMax.X) then Swap(KMin, KMax);
  if (KMid.X > KMax.X) then Swap(KMid, KMax);
  if ((KMin.X = KMid.X) and (KMin.Y > KMid.Y)) then Swap(KMin, KMid);
  if ((KMid.X = KMax.X) and (KMid.Y > KMax.Y)) then Swap(KMid, KMax);
  
  //KNew := new Point3D(KMid.X, (KMax.Y - KMin.Y) / (KMax.X - KMin.X) * (KMid.X - KMax.X) + KMax.Y, 0);
  KNew := new Point3D(KMid.X, (KMax.Y - KMin.Y) / (KMax.X - KMin.X) * (KMid.X - KMin.X) + KMin.Y, 0);
  //KNew := new Point3D(KMid.X, (KMax.Y-KMin.Y)/(KMax.X - KMin.X)*KMid.X + KMin.Y - (KMax.Y-KMin.Y)/(KMax.X - KMin.X)*KMin.X, 0); // Это одно и то же
  
  // Левые пределы для интегрирования по треугольнику
  FxLU := Fx(KMin, KNew);
  FxLD := Fx(KMin, KMid);
  if (KMid.Y > KNew.Y) then Swap(FxLU, FxLD);
  
  // Правые пределы для интегрирования по треугольнику
  FxRU := Fx(KNew, KMax);
  FxRD := Fx(KMid, KMax);
  if (KMid.Y > KNew.Y) then Swap(FxRU, FxRD);
  
  
  // Коэффиценты при базисе (K0-K1) вида (a0x+b0y+c0)/d0
  d0 := (K2.X - K1.X) * (K1.Y - K2.Y) - (K0.X - K1.X) * (K1.Y - K2.Y) - (K0.Y - K2.Y) * (K2.X - K1.X); 
  a0 := K2.Y - K1.Y;
  b0 := K1.X - K2.X;
  c0 := (K2.X - K1.X) * (K1.Y - K2.Y) + K1.X * (K1.Y - K2.Y) + K2.Y * (K2.X - K1.X);
  
  // Коэффиценты при базисе (K1-K2) вида (a1x+b1y+c1)/d1
  d1 := (K0.X - K2.X) * (K2.Y - K0.Y) - (K1.X - K2.X) * (K2.Y - K0.Y) - (K1.Y - K0.Y) * (K0.X - K2.X); 
  a1 := K0.Y - K2.Y;
  b1 := K2.X - K0.X;
  c1 := (K0.X - K2.X) * (K2.Y - K0.Y) + K2.X * (K2.Y - K0.Y) + K0.Y * (K0.X - K2.X);
  
  // Коэффиценты при базисе (K2-K0) вида (a2x+b2y+c2)/d2
  d2 := (K1.X - K0.X) * (K0.Y - K1.Y) - (K2.X - K0.X) * (K0.Y - K1.Y) - (K2.Y - K1.Y) * (K1.X - K0.X); 
  a2 := K1.Y - K0.Y;
  b2 := K0.X - K1.X;
  c2 := (K1.X - K0.X) * (K0.Y - K1.Y) + K0.X * (K0.Y - K1.Y) + K1.Y * (K1.X - K0.X);
  
  
  var EEE := 1.0E10;
  if ((Abs(d0) > EEE) or (Abs(d1) > EEE) or (Abs(d2) > EEE)) then Println(T0.N, ':', d0, d1, d2);
  
  var EEE_ := 1.0E-10;
  if ((Abs(d0) < EEE_) or (Abs(d1) < EEE_) or (Abs(d2) < EEE_)) then Println(T0.N, ':', d0, d1, d2);
  
  
  N0 := XYtoN(K0) - r_dt_nt[1][1];
  N1 := XYtoN(K1) - r_dt_nt[1][1];
  N2 := XYtoN(K2) - r_dt_nt[1][1];
 
  //Println(N0,N1,N2);
  
  // N0
  if N0 >= 0 then begin
    u_arr.Add_U(integrate_basis(a0, a0, b0, b0, c0, c0, d0, d0, true),N0,N0);
    b_arr[N0] += integrate_basis(a0, a0, b0, b0, c0, c0, d0, d0, false) * T0.P0.Z;
  end;
  
  if (N0 >= 0) and (N1 >= 0) then begin
    u_arr.Add_U(integrate_basis(a0, a1, b0, b1, c0, c1, d0, d1, true),N0,N1);
    b_arr[N0] += integrate_basis(a0, a1, b0, b1, c0, c1, d0, d1, false) * T0.P1.Z;
  end;
  
  if (N0 >= 0) and (N2 >= 0) then begin
    u_arr.Add_U(integrate_basis(a0, a2, b0, b2, c0, c2, d0, d2, true),N0,N2);
    b_arr[N0] += integrate_basis(a0, a2, b0, b2, c0, c2, d0, d2, false) * T0.P2.Z;
  end;
  
  // N1
  
  if (N1 >= 0) and (N0 >= 0) then begin
    u_arr.Add_U(integrate_basis(a1, a0, b1, b0, c1, c0, d1, d0, true),N1,N0);
    b_arr[N1] += integrate_basis(a1, a0, b1, b0, c1, c0, d1, d0, false) * T0.P0.Z;
  end;
  
  if (N1 >= 0) then begin
    u_arr.Add_U(integrate_basis(a1, a1, b1, b1, c1, c1, d1, d1, true),N1,N1);
    b_arr[N1] += integrate_basis(a1, a1, b1, b1, c1, c1, d1, d1, false) * T0.P1.Z;
  end;
  
  if (N2 >= 0) and (N1 >= 0) then begin
    u_arr.Add_U(integrate_basis(a1, a2, b1, b2, c1, c2, d1, d2, true),N1,N2);
    b_arr[N1] += integrate_basis(a1, a2, b1, b2, c1, c2, d1, d2, false) * T0.P2.Z;
  end;
  
  
  // N2
  if (N2 >= 0) and (N0 >= 0) then begin
    u_arr.Add_U(integrate_basis(a2, a0, b2, b0, c2, c0, d2, d0, true),N2,N0);
    b_arr[N2] += integrate_basis(a2, a0, b2, b0, c2, c0, d2, d0, false) * T0.P0.Z;
  end;
  
  if (N2 >= 0) and (N1 >= 0) then begin
    u_arr.Add_U(integrate_basis(a2, a1, b2, b1, c2, c1, d2, d1, true),N2,N1);
    b_arr[N2] += integrate_basis(a2, a1, b2, b1, c2, c1, d2, d1, false) * T0.P1.Z;
  end;
  
  if (N2 >= 0) then begin
    u_arr.Add_U(integrate_basis(a2, a2, b2, b2, c2, c2, d2, d2, true),N2,N2);
    b_arr[N2] += integrate_basis(a2, a2, b2, b2, c2, c2, d2, d2, false) * T0.P2.Z;
  end;
  
end;

procedure SetCamOrt();
begin
 Graph3D.hvp.Camera := new System.Windows.Media.Media3D.OrthographicCamera;
end;

Begin
  SetConsoleIO;
  GraphABC.Window.Minimize;
  Graph3D.View3D.ShowGridLines := false;
  Lights.AddPointLight(GrayColor(55), new Point3D(50000, -50000, -100000));
  
  View3D.SubTitle := 'F - F(x)'#10'U - U(x)'#10'D - триангуляция Делоне'#10'N - нумерация'#10'T - медленная тестовая нумерация'#10'B - базисы';
  
  //Graph3D.View3D.HideAll;
  //Graph3D.ShowCoordinateSystem := false;
  
  Nsum := Delone(R, dr);
  
  Println('dn', R/dr);
  Println('Nsum', Nsum, XYtoN(new Point3D(0.0, 0.0, 0.0)));
  Println('Вершины первого ряда:', r_dt_nt[1][1]);
  Println('Размер матрицы:', Nsum-r_dt_nt[1][1]);
  Println('Всего вершин:', Nsum);
  
  Nsum -= r_dt_nt[1][1];
  SetLength(b_arr,Nsum);
  u_arr := new M_Zip;
  
  //ZfromPic('t1.jpg'); // leon.jpg skull.jpg  t1.jpg
  ZfromFunc;
  
  var DT := DateTime.Now;
  
  Println('Start integrate');
  DT := DateTime.Now;
  for var k := 0 to TrArr.Count - 1 do
    for var j := 0 to TrArr[k].Count - 1 do
      integral_triangle(TrArr[k][j]);
  Println('End integrate:', DateTime.Now - DT);
  
  //Print(u_arr.M_rows);
  
  DT := DateTime.Now;
  Println('Start solving');
  u_arr.Solve(b_arr);
  //var c_arr := u_arr.Inv * b_arr;
  Println('End solving:', DateTime.Now - DT);
  var c_arr := u_arr.X_vec;
  
  DT := DateTime.Now;
  Println('Start equlaizing');
  DT := DateTime.Now;
  for var k := 0 to TrArr.Count - 1 do
    for var j := 0 to TrArr[k].Count - 1 do
    begin
      try
        TrArr[k][j].U0.Z := c_arr[XYtoN(TrArr[k][j].U0) - r_dt_nt[1][1]];
      except end;
      try
        TrArr[k][j].U1.Z := c_arr[XYtoN(TrArr[k][j].U1) - r_dt_nt[1][1]];
      except end;
      try
        TrArr[k][j].U2.Z := c_arr[XYtoN(TrArr[k][j].U2) - r_dt_nt[1][1]];
      except end;
      
      // Для расчета Оси Z на графике
      var Zmin, Zmax: real;
      Zmin := Min(TrArr[k][j].U0.Z, TrArr[k][j].U1.Z, TrArr[k][j].U2.Z);
      Zmax := Max(TrArr[k][j].U0.Z, TrArr[k][j].U1.Z, TrArr[k][j].U2.Z);
      if Z_Ax_min > Zmin then Z_Ax_min := Zmin;
      if Z_Ax_max < Zmax then Z_Ax_max := Zmax;
      
      TrArr[k][j].D0.Z := TrArr[k][j].P0.Z - TrArr[k][j].U0.Z;
      TrArr[k][j].D1.Z := TrArr[k][j].P1.Z - TrArr[k][j].U0.Z;
      TrArr[k][j].D2.Z := TrArr[k][j].P2.Z - TrArr[k][j].U0.Z;
    end;
  Println('End equlaizing:', DateTime.Now - DT);
  
  
  DT := DateTime.Now;
  Println('Start draw');
  
  Graph3D.OnKeyDown := KK -> begin
    if IsDraw then 
      exit;
    Graph3D.Window.Title := 'Формируем рисунок';
    case KK of
      Key.U: Redraw(Draw3D_U);
      Key.F: Redraw(Draw3D_F);
      Key.D: Redraw(Draw3D_D); // Делоне
      Key.N: Redraw(Draw3D_N); // Нумерация
      Key.T: Redraw(Draw3D_TN); // Медленная тестовая нумерация
      Key.B: Redraw(Draw3D_B); // Базисы в "0" треугольнике
    end;
    Graph3D.Window.Title := 'Рисунок сформирован';
  end;
  

  
  // Позиционируем камеру
  Invoke(SetCamOrt);
  //Graph3D.View3D.ShowCameraInfo := true;
  Graph3D.View3D.CameraMode := CameraMode.Inspect;
  //Camera.UpDirection := V3D(0, 0, 1);
  //Camera.Rotate(new Vector3D(1, 0, 0),10);
  Camera.Rotate(V3D(0, 1, 0),-45);
  Camera.UpDirection := V3D(0,0,1);
  Camera.Position := P3D(5, 5, 0);
  Camera.LookDirection := V3D(-5, -5, 0);
  

//  Camera.Distanse := 0.5; // Расстояние до "среза" рендеренга
//  Camera.UpDirection := V3D(0,0,1);
  //Graph3D.View3D.CameraMode := CameraMode.WalkAround;
  
  Graph3D.Window.Title := 'Формируем рисунок';
  // Рисуем Оси координат
  Redraw(Draw3D_Ax);
  
  // Рисуем U(x)
  Redraw(Draw3D_U);
  
  
  Graph3D.Window.Title := 'Рисунок сформирован';
  Println('End draw:', DateTime.Now - DT);
  DT := DateTime.Now;
  
  
  var dotnow: Point3D;
  Graph3D.OnMouseMove += procedure (x, y, mb) -> begin
    dotnow := FindNearestObjectPoint(x, y);
    if dotnow = BadPoint then begin
      Graph3D.Window.Caption := '--- --- ---';
    end
    else begin
      Graph3D.Window.Caption := ((dotnow.X*1000).Round / 1000).ToString + ' ' + ((dotnow.Y*1000).Round / 1000).ToString + ' ' + ((dotnow.Z*1000).Round / 1000).ToString;
    end;
  end;
  
  Println('U(0,0,0) =', c_arr.Last);
  
  //u_arr.Println;
  // Расхождение СЛАУ
//  Println;
//  Println('Строки СЛАУ с расхождением |M(u)- b| > 1.0E-10 :');
//  var err: real;
//  for var i := 1 to u_arr.Count do
//  begin
//    err := u_arr.Row(i) * c_arr - b_arr.Value[i - 1];
//    if Abs(err) > 1.0E-10 then Println(i, err);
//  end;
  Println('THE END.');
//  // Чистим память
//  b_arr.Length := 0;
//  u_arr.ColCount := 0;
//  u_arr.RowCount := 0;
end.