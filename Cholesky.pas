unit Cholesky;

//uses NumLibABC;

type
  a_Zip = class(IComparable<a_Zip>)
  private
  
  public
    U: real;
    L: real;
    i: integer;
    constructor(i_col: integer; U_value: real;  L_value: real := 0);
    begin
      i := i_col;
      U := U_value;
      L := L_value;
    end;
    
    function CompareTo(b: a_Zip): integer;
    begin
      Result := i - b.i;
    end;
  end;

type
  i_Zip = List<a_Zip>;

type
  M_Zip = class
  private
    NN: integer := 0;
    Lt := new List<List<(integer, real)>>; // Траспонированная "сжатая" матрица
    Y_vec := new List<real>;
    M_rows := new List<i_Zip>; // Основная матрица
    
  public
    X_vec := new List<real>;
    
    function Count() := NN;
    
    procedure Add_U(U: real; j, i: integer);
    var
      count := 0;
    begin
    
      if j < 0 then exit;
      if i < 0 then exit;
      if j < i then exit;
      if U = 0 then exit;
      
      
      for var ii := 0 to j-M_rows.Count do M_rows.Add(new i_Zip);
      if M_rows[j].Count = 0 then begin M_rows[j].Add(new a_Zip(i, U)); NN := M_rows.Count; exit; end;
      
      foreach var S in M_rows[j] do
        if S.i < i then Inc(count)
        else break;
      
      if ((count >= M_rows[j].Count) or (M_rows[j][count].i <> i)) then M_rows[j].Insert(count, new a_Zip(i, U))
      else begin 
        M_rows[j][count].U += U; 
        if M_rows[j][count].U = 0 then M_rows[j].RemoveAt(count);
      end;
      
      NN := M_rows.Count;
    end; 
    
    
    procedure Add_L(L: real; j, i: integer);
    var
      count := 0;
    begin
      
      if j < 0 then exit;
      if i < 0 then exit;      
      if j < i then exit;
      if L = 0 then exit;
      
      foreach var S in M_rows[j] do
        if S.i < i then Inc(count)
        else break;
      
      if ((count >= M_rows[j].Count) or (M_rows[j][count].i <> i)) then M_rows[j].Insert(count, new a_Zip(i, 0, L))
      else begin
        M_rows[j][count].L += L;
        if (M_rows[j][count].U = 0) and (M_rows[j][count].L = 0) then M_rows[j].RemoveAt(count);
      end;
    end; 
    
    
    /// Разбиение матрицы LU методом Холецкого
    procedure LLt();
    begin
      var sum := 0.0;
      var ki := 0;
      for var i := 0 to NN do
      begin
        for var j := 0 to i - 1 do
        begin
          sum := 0;
          ki := 0;
          for var k := 0 to M_rows[j].Count - 2 do
          begin // Накапливаем сумму. Не используем Li и Ui т.к. они медленные - каждый раз перебирают всю строку сначала
            while M_rows[i][ki].i < M_rows[j][k].i do Inc(ki); // Ищем следующий U(i,k) у которого .i совпадает с U(j,k)
            if M_rows[i][ki].i = M_rows[j][k].i then
              sum += M_rows[i][ki].L * M_rows[j][k].L;
          end;
          var Uij := 0.0;
          while M_rows[i][ki].i < M_rows[j][M_rows[j].Count - 1].i do Inc(ki); // "Дотягиваем" счетчик ki до j
          if M_rows[i][ki].i = M_rows[j][M_rows[j].Count - 1].i then Uij := M_rows[i][ki].U; // U(i,ki) не "0" только если имеет тот же .i, что U(j,j)
          Add_L((Uij - sum) / M_rows[j][M_rows[j].Count - 1].L, i, j); // U(j,j) последний в своей строке
        end;
        
        // Экономим на постоянной проверке в вышестоящем цикле, т.к. понятно, что i=j только на последней итерации
        sum := 0;
        ki := 0;
        for var k := 0 to M_rows[i].Count - 2 do
        begin // Аналогично вышеизложенному для суммы
          while (M_rows[i][ki].i < M_rows[i][k].i) do Inc(ki);
          if M_rows[i][ki].i = M_rows[i][k].i then 
            sum += M_rows[i][ki].L * M_rows[i][k].L;
        end;
        Add_L(Sqrt(M_rows[i][M_rows[i].Count - 1].U - sum), i, i); // Не используем Ui т.к. ежу понятно, что U(i,i) последний член в этом ряду
      end;
    end;
    
    /// Решение для нижней треугольной матрицы Ly=B => Y_vec + формируем Lt
    procedure L_solve(Bb: array of real);
    begin
      Y_vec.Add(Bb[0] / M_rows[0][0].L);
      
      for var j := 1 to NN do
      begin
        var sum := 0.0;
        for var i := 0 to M_rows[j].Count - 2 do 
          sum += M_rows[j][i].L * Y_vec[M_rows[j][i].i]; // Накапливаем сумму для вычисления Yj
        Y_vec.Add((Bb[j] - sum) / M_rows[j][M_rows[j].Count - 1].L);
      end;
    end;
    
    
    /// Транспонируем L в Lt в "сжатом" формате
    procedure Lt_make();
    begin
      // Добавляем в Lt NN строк, т.к. все они точно нужны
      loop M_rows.Count do
        Lt.Add(new List<(integer, real)>);
      // Для решения СЛАУ нам не важен порядок хранения недиагональных элементов в строке Lt, т.к. свой i они хранят в себе,
      // а диагональные понятно как "положить" в строку последними, сделаем это позже
      for var j := 1 to NN do // Обходим без диагональных!!!
        for var i := 0 to M_rows[j].Count - 2 do
          Lt[M_rows[j][i].i].Add((j, M_rows[j][i].L));
        
      // Добавляем диагональные элементы в конец - только для компактности записи решения СЛАУ, 
      // это лишняя трата времени и памяти
      for var j := 0 to NN do
        Lt[j].Add((j, M_rows[j][M_rows[j].Count - 1].L));
    end;
      
      
    /// Решение для верхней треугольной матрицы Ux=Y => X_vec
    procedure Lt_solve();
    begin
      Lt_make(); // Транспонируем матрицу L
      
      // Решаем СЛАУ обратным ходом снизу-вверх
      X_vec.Add(Y_vec[NN] / Lt[NN][0].Item2);
            for var j := NN-1 downto 0 do
            begin
              var sum := 0.0;
              for var i := 0 to Lt[j].Count - 2 do
                sum += Lt[j][i].Item2 * X_vec[X_vec.Count-1 - NN + Lt[j][i].Item1]; // Накапливаем сумму для вычисления Yj
              X_vec.Insert(0,(Y_vec[j] - sum) / Lt[j][Lt[j].Count - 1].Item2);
            end;
    end;
    
    procedure PrintM_U(w: integer := 6; d: integer := 2; fl: boolean := true);
    begin
      Println;
      Println('Матрица', NN+1, 'x', NN+1,'по U:');
      for var i := 0 to NN do write(1.0*i:w:0);
      Println;
      for var j := 0 to NN do begin
        var ki := 0;
        for var i := 0 to NN do
          if M_rows[j][ki].i = i then begin 
            write(M_rows[j][ki].U:w:d); //M_rows[j][ki].U
            if ki < M_rows[j].Count-1 then Inc(ki);
          end
          else write(0.0:w:d);
        if fl then
          Println ('  >', j)
        else Println();
      end;
      Println;
    end;
    
    procedure PrintM_L(w: integer := 6; d: integer := 2; fl: boolean := true);
    begin
      Println;
      Println('Матрица', NN+1, 'x', NN+1,'по L:');
      for var i := 0 to NN do write(1.0*i:w:0);
      Println;
      for var j := 0 to NN do begin
        var ki := 0;
        for var i := 0 to NN do
          if (M_rows[j][ki].i = i) then begin 
            write(M_rows[j][ki].L:w:d); //M_rows[j][ki].U
            if ki < M_rows[j].Count-1 then Inc(ki);
          end
          else write(0.0:w:d);
        if fl then
          Println ('  >', j)
        else Println();
      end;
      Println;
    end;
    
    procedure PrintRowsM(fl:boolean := true);
    begin
      Println;
      Println('Матрица', NN+1, 'x', NN+1,'по U:');
      Println;
      for var j := 0 to NN do begin
        M_rows[j].Print;
        if fl then
          Println ('  >', j)
        else Println();
      end;
    end;
    
//    function UnZip_U(): Matrix;
//    begin
//    var Mtemp := new Matrix(NN+1,NN+1);
//      for var j := 0 to NN do
//        for var i := 0 to M_rows[j].Count-1 do begin
//          Mtemp.Value[j,M_rows[j][i].i] := M_rows[j][i].U;
//          Mtemp.Value[M_rows[j][i].i,j] := Mtemp.Value[j,M_rows[j][i].i];
//        end;
//      Result := Mtemp;
//    end;
//    
//    function UnZip_L(): Matrix;
//    begin
//    var Mtemp := new Matrix(NN+1,NN+1);
//      for var j := 0 to NN do
//        for var i := 0 to M_rows[j].Count-1 do
//          Mtemp.Value[j,M_rows[j][i].i] := M_rows[j][i].L;
//      Result := Mtemp;
//    end;
//    
//    function UnZip_Lt(): Matrix;
//    begin
//    var Mtemp := new Matrix(NN+1,NN+1);
//      for var j := 0 to NN do
//        for var i := 0 to Lt[j].Count-1 do 
//          Mtemp.Value[j,Lt[j][i].Item1] := Lt[j][i].Item2;
//      Result := Mtemp;
//    end;
    
    
        
    procedure Solve(Bb: array of real);
    begin
      var DT := DateTime.Now;
      var ttt := new Stopwatch();
      NN := M_rows.Count-1;
      
      ttt.Start;
      LLt();
      Println('End LLt:', DateTime.Now - DT, ttt.ElapsedMilliseconds, ttt.ElapsedTicks);
      
      DT := DateTime.Now;
      ttt.Reset;
      ttt.Start;
      L_solve(Bb);
      Println('End L_solve:', DateTime.Now - DT, ttt.ElapsedMilliseconds, ttt.ElapsedTicks);
      
      DT := DateTime.Now;
      ttt.Reset;
      ttt.Start;
      Lt_solve();
      Println('End Lt_solve:', DateTime.Now - DT, ttt.ElapsedMilliseconds, ttt.ElapsedTicks);
    end;
  end;
end.