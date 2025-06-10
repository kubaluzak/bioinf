# Sekwencjonowanie DNA przez hybrydyzację w oparciu o chipy binarne (z uwzględnieniem błędów pozytywnych)

---

## Autorzy  
**Jakub Urbaniak 155870**  
**Wiktor Wachowski 155859**

---

## Spis treści

0. [Wprowadzenie i opis problemu](#0-opis-i-formalizacja-problemu)  
1. [Wprowadzenie i opis problemu](#1-wprowadzenie-i-opis-problemu)  
2. [Algorytm BFS — dokładny](#2-algorytm-bfs--dokładny)  
   2.1. [Reprezentacja stanu](#21-reprezentacja-stanu)  
   2.2. [Inicjalizacja BFS](#22-inicjalizacja-bfs)  
   2.3. [Iteracyjny proces BFS](#23-iteracyjny-proces-bfs)  
   2.4. [Zakończenie i wynik](#24-zakończenie-i-wynik)  
3. [Analiza złożoności czasowej](#3-analiza-złożoności-czasowej)  
   3.1. [Parametry](#31-parametry)  
   3.2. [Przestrzeń stanów](#32-przestrzeń-stanów)  
   3.3. [Koszt pojedynczego stanu](#33-koszt-pojedynczego-stanu)  
   3.4. [Całkowita złożoność](#34-całkowita-złożoność)  
4. [Zależność czasu działania algorytmu od długości sekwencji](#4-zależność-czasu-działania-algorytmu-od-długości-sekwencji)

---

## 0. Wprowadzenie i opis problemu

W niniejszym sprawozdaniu pochylimy się nad problemem sekwencjonowania DNA, czyli odtwarzania pełnej sekwencji nukleotydów organizmu na podstawie danych eksperymentalnych. Jednym z podejść do tego zagadnienia jest sekwencjonowanie przez hybrydyzację, które polega na wykorzystywaniu krótkich fragmentów DNA, zwanych sondami, do odtwarzania oryginalnego ciągu.

W naszym projekcie zajmujemy się wariantem tego problemu z wykorzystaniem chipów binarnych. W tej metodzie każda sonda nie wskazuje bezpośrednio konkretnego fragmentu DNA, lecz jest zakodowana w postaci binarnej, co wprowadza niejednoznaczność i utrudnia odtworzenie oryginalnej sekwencji.
Dane dotyczące sond reprezentowane są przez spektra, które kategoryzują nukleotydy opierając się na ich różnych właściwościach.

Dodatkową trudnością w tym wariancie są błędy pozytywne. Są to nadmiarowe sondy, które nie znajdowały się w pierwotnej sekwencji DNA. Ich obecność jest efektem niedoskonałości technologicznych lub zakłóceń w procesie hybrydyzacji. W praktyce oznacza to, że nie każda dopasowana sonda jest poprawnym fragmentem DNA.

Ostatecznie, biorąc pod uwagę te aspekty, zadaniem algorytmu jest odtworzenie pełnej sekwencji DNA o zadanej długości, rozpoczynając od znanego fragmentu startowego. Proces budowania sekwencji odbywa się nukleotyd po nukleotydzie - na każdym kroku wyszukiwane są sondy, które pasują do aktualnie utworzonego fragmentu w obu spektrach, a następnie wybierane jest takie rozszerzenie, które zgodnie łączy oba kodowania i najlepiej spełnia warunki spójności. Dopiero złożenie dwóch sond z przeciwnych spektrów zawęża możliwe wartości do konkretnego nukleotydu, umożliwiając stopniowe odtwarzanie oryginalnego ciągu.

---

## 1. Formalizacja problemu

### Dane wejściowe
- **Długość badanego fragmentu DNA (n)**, do którego chcemy zrekonstruować pełną sekwencję
- **Długość sond oligonukleotydowych (k)**
- **Sekwencja startowa (S<sub>0</sub>)** o długości *k* — znany fragment DNA, od którego rozpoczynamy rekonstrukcję całej sekwencji,  
- **Sondy** w dwóch spektrach:  
   - **Spektrum WS (Weak/Strong)** - dzieli nukleotydy według siły ich wiązań wodorowych: **A,T -> W** (wiązanie podwójne); **C,G -> S** (wiązanie potrójne)
   - **Spektrum RY (Purine/Pyrimidine)** - dzieli nukleotydy według ich typu zasady azotowej: **A,G -> R** (puryny); **C,T -> Y** (pirymidyny)

> [!NOTE]
> *Ostateczna długość sond wynosi w praktyce **2k−1**, ponieważ każda pozycja w sekwencji - poza ostatnią - jest reprezentowana przez spektra WS i RY. Sonda z danego spektrum nosi ze sobą połowę informacji o nukleotydzie. Dopiero analiza wspólna WS i RY pozwala jednoznacznie określić konkretny nukleotyd, dlatego w praktyce pozycji w sondzie jest więcej.*

### Cel

Celem problemu jest znalezienie sekwencji *R*, która:
- Odtwarza pełną sekwencję DNA o zadanej długości n
- Jest superciągiem względem zbioru sond — tzn. fragmenty sond zachodzą na siebie na k-1 pozycjach (suffix poprzedniego wyrazu = prefix następnego wyrazu).
  Formalnie:
           dla każdej pary kolejnych sond *S<sub>y</sub>* i *S<sub>x</sub>*, gdzie *S<sub>y</sub>* występuje przed *S<sub>x</sub>*, zachodzi:
                                                         *n<sub>y+i</sub>* = *n<sub>x+i-1</sub>*
- Rozszerza znany fragment startowy *S<sub>0</sub>*,  
- Jest spójny z dostępnymi sondami w obu spektrach (dla każdego dopisanego nukleotydu istnieją sondy w spektrum WS i RY o odpowiednich prefixach i dopisanych nukleotydach),  
- Minimalizuje błędy pozytywne wynikające z obecnośći nadmiarowych sond,

---

## 2. Części wspólne algorytmów

Zarówno algorytm dokładny, jak i heurystyczny bazują na wspólnych mechanizmach przetwarzania danych wejściowych oraz przygotowania struktur pomocniczych wykorzystywanych do dalszego działania algorytmów. Poniżej przedstawiono elementy wspólne dla obu implementacji:

### 2.1. Wczytywanie danych wejściowych z plików XML

W obu przypadkach algorytm rozpoczyna swoje działanie od wczytania danych z pliku XML. Każdy plik zawiera:
- startową sekwencję
- docelową długość sekwencji
- dwa zbiory sond (dla każdego ze spektrów)

### 2.2. Wczytywanie danych wejściowych z plików XML

Oba algorytmy korzystają z metod kodowania sekwencji na dwa sposoby:
WS (A i T → W; C i G → S)
RY (A i G → R; C i T → Y)

### 2.3. Budowanie map prefiksów sond

W celu efektywnego wyszukiwania możliwych dalszych nukleotydów do budowania potencjalnej sekwencji, w obu algorytmach tworzona jest mapa prefiksów. Kluczami map są prefiksy długości *k-1* danej sondy (bez ostatniego znaku), a wartościami są zbiory liter (nukleotydów) możliwych do dopisania jako kolejnych.

---

## 3. Algorytm BFS — dokładny

Do rozwiązania problemu wykorzystujemy algorytm **przeszukiwania wszerz (BFS)** po przestrzeni stanów reprezentującej możliwe rozszerzenia rekonstrukcji sekwencji.

### 3.1. Reprezentacja stanu

Każdy stan algorytmu reprezentowany jest przez:

- **Suffix sekwencji** z ostatnich \(k-1\) nukleotydów zakodowanych w spektrum WS,  
- **Suffix sekwencji** z ostatnich \(k-1\) nukleotydów zakodowanych w spektrum RY,  
- Aktualnie zrekonstruowaną sekwencję nukleotydów (ciąg liter A, C, T, G).

W kodzie, zamiast przechowywać całe sekwencje, operujemy na suffixach i rozszerzamy je o kolejne nukleotydy zgodne z sondami.

### 3.2. Inicjalizacja BFS

- Na początku kodujemy suffix startowej sekwencji w obu spektrach,  
- Tworzymy dwie mapy prefixów (po jednym na każde spektrum), które dla prefixu długości \(k-1\) podają zbiór możliwych nukleotydów do dopisania,  
- Inicjujemy kolejkę BFS z początkowym stanem: suffixy startowe i startowa sekwencja,  
- Tworzymy zbiór `visited` by pamiętać, które pary suffixów już zostały przetworzone — to zapobiega powtórkom.

### 3.3. Iteracyjny proces BFS

W każdej iteracji:

1. Pobieramy stan z kolejki — aktualne suffixy WS i RY oraz dotychczas zrekonstruowaną sekwencję.  
2. Sprawdzamy, czy długość sekwencji osiągnęła wymaganą wartość — jeśli tak, to zwracamy rozwiązanie.  
3. Pobieramy z map prefixów kandydujące nukleotydy, które można dołączyć do aktualnych suffixów WS i RY.  
4. Iterujemy przez nukleotydy \( \{A, C, T, G\} \) i sprawdzamy, czy nukleotyd pasuje jednocześnie do obu spektrów (czy występuje jako dozwolony w obu mapach dla odpowiednich prefixów).  
5. Dla każdego pasującego nukleotydu tworzymy nowy stan:  
   - Obliczamy nowe suffixy, przesuwając istniejący suffix o 1 pozycję i dopisując kodowaną literę (w spektrum WS i RY),  
   - Doklejamy nukleotyd do aktualnej sekwencji,  
   - Jeśli para suffixów nie była wcześniej odwiedzona, dodajemy nowy stan do kolejki i do zbioru `visited`.

### 3.4. Zakończenie i wynik

- Algorytm działa dopóki kolejka nie jest pusta lub nie znajdzie sekwencji o wymaganej długości,  
- W przypadku braku rozwiązania informuje o tym użytkownika.

---

## 3. Analiza złożoności czasowej

### 3.1. Parametry

| Symbol | Znaczenie                       |
|--------|--------------------------------|
| \(k\)  | Długość sondy                  |
| \(n\)  | Długość docelowej sekwencji   |
| \(S\)  | Liczba sond w każdym spektrum |

### 3.2. Przestrzeń stanów

Algorytm operuje na parach suffixów długości \(k-1\) z obu spektrów.

Liczba stanów do rozpatrzenia:

\[
O(S^2)
\]

### 3.3. Koszt pojedynczego stanu

- Pobranie kandydatów i sprawdzenie odwiedzenia to operacje w czasie stałym \(O(1)\).

### 3.4. Całkowita złożoność

\[
O(S^2)
\]

z założeniem, że eksplorowane są wszystkie możliwe pary prefixów. W praktyce jest zwykle mniejsza dzięki filtrowaniu i zatrzymaniu po znalezieniu rozwiązania.

---

## 4. Algorytm Beam Search - heurystyczny

## Opis algorytmu

Algorytm Beam Search jest heurystyczną modyfikacją algorytmu BFS, stosowaną w celu ograniczenia eksploracji przestrzeni rozwiązań. Kluczowym założeniem algorytmu jest to, że w każdej iteracji rozwijanych jest jedynie maksymalnie *X* potencjalnie najlepszych kandydatów. Kandydaci wybierani są na podstawie priorytetu, który jest wyliczany heurystycznie - uwzględnia długość obecnej sekwencji oraz częstotliwość występowania nukleotydów w dostępnych sondach. Dzięki temu liczba rozważanych ścieżek pozostaje kontrolowana, co pozwala na szybsze działanie kosztem potencjalnego pominięcia najlepszego rozwiązania.

### 4.1. Reprezentacja stanu

Analogicznie jak w algorytmie dokładnym (BFS), każdy stan algorytmu reprezentowany jest przez:

- **Suffix sekwencji** z ostatnich *k-1* nukleotydów zakodowanych w spektrum WS,  
- **Suffix sekwencji** z ostatnich *k-1* nukleotydów zakodowanych w spektrum RY,  
- Aktualnie zrekonstruowaną sekwencję nukleotydów (ciąg liter alfabetu A, C, T, G).

Stany przechowywane są natomiast w kopcu priorytetowym, gdzie priorytet określany jest heurystycznie na podstawie częstotliwości danego nukleotydu w zbiorze sond oraz długości sekwencji.

### 4.2. Kroki działania algorytmu

Algorytm składa się z następujących faz:

1. Odczytujemy z pliku parametry wejściowe. Dla obu zbiorów sond budowane są mapy prefiksów, ułatwiające wyszukiwanie możliwych rozszerzeń sekwencji (podobnie jak w BFS).
2. Inicjalizacja kolejki priorytetowej (kopca) - algorytm umieszcza w kolejce startowy stan, obejmujący zakodowane suffiksy WS i RY oraz startową sekwencję. Każdemu stanowi przypisywany jest priorytet wyliczany na podstawie długości sekwencji oraz częstości liter w sondach.
3. W części iteracyjnej programu:
   - W każdej iteracji rozwijamy maks *X* najlepszych (z najwyższym priorytetem) kandydatów z kolejki
   - Dla każdego rozwijanego stanu algorytm generuje nowe kandydatury, dodając kolejne nukleotydy, pod warunkiem że są zgodne z sondami WS i RY (na podstawie map prefiksów)
   - Każdy nowy stan dodawany jest do kolejki priorytetowej jeśli nie został jeszcze odwiedzony
4. Jeśli podczas rozwijania stanów uda się uzyskać sekwencję o żądanej długości, algorytm zwraca wynik. W przeciwnym wypadku, gdy kolejka się wyczerpie, algorytm kończy działanie bez znalezienia rozwiązania.

#### Kroki szczegółowe dla poszczególnego stanu
1. Pobieramy listę możliwych nukleotydów do dopisania dla aktualnego suffiksu w obu spektrach
2. Wyznaczamy część wspólną obu zbiorów kandydatów - tylko nukleotydy akceptowalne przez oba spektra
3. Jeśli istnieje tylko jeden kandydat, dopisujemy go do sekwencji
4. Jeśli istnieje więcej niż jeden kandydat, dla każdego z nich tworzony jest nowy stan z wyliczonym priorytetem na podstawie długości sekwencji oraz zgodności z sondami. Priorytet ten określa kolejność rozwijania stanów w kolejnych iteracjach.
5. Jeśli natomiast nie ma żadnego możliwego dopisania, kończymy rozwijanie tego stanu (bez dodawania nowego do kolejki)
6. Aktualizujemy suffix i idziemy do kolejnego kroku

### 4.3. Analiza złożoności algorytmu Beam Search

**Parametry:**
- **`k`** — długość sondy (długość suffixu to `k−1`),  
- **`n`** — docelowa długość rekonstruowanej sekwencji,  
- **`S`** — liczba sond w każdym spektrum,  
- **`X`** — szerokość wiązki (beam width), czyli maksymalna liczba kandydatów rozwijanych w każdej iteracji  

**Analiza:**
1. W każdej iteracji rozwijanych jest maksymalnie **`X`** najlepszych stanów (kandydatów)  
2. Dla każdego stanu generujemy do **4 nowych kandydatów** (A, C, G, T), z których dopuszczamy tylko te zgodne z oboma spektrami  
3. Liczba iteracji wynosi maksymalnie **`n − |start|`**, gdzie `|start|` to długość sekwencji startowej.  

**Złożoność czasowa:**  
`O(X ⋅ n)`  
(co oznacza liniowy czas względem długości docelowej sekwencji i szerokości wiązki).  

**Złożoność pamięciowa:**  
`O(X ⋅ n)`  
(z uwagi na przechowywanie do `X` sekwencji o długości `n` oraz ich suffixów).  

> [NOTE]
> W zapisie złożoności czasowej i pamięciowej używamy proszczonego symbolu *n*, mimo że w rzeczywistości liczba maks. iteracji algorytmów jest mniejsza = *n - |start|*. Wynika to z tego, że różnica jest zaniedbywalna, zwłaszcza dla większych danych testowych.

---

## 5 Analiza i porównanie wyników

## 5.1. Ograniczenia w danych testowych

Należy wpierw zwrócić uwagę na pewne niedogodnienia związane z procesem testowania algorytmów:

Głównym ograniczeniem przeprowadzonych testów był stosunkowo wąski zakres parametrów wejściowych, który nie pozwolił w pełni uwidocznić różnic w działaniu obu algorytmów. Ze względów technicznych eksperymenty ograniczono do sekwencji o długości nieprzekraczającej 1000 nukleotydów, choć teoretyczne możliwości implementacyjne pozwalały na obsługę sekwencji do 65,535 elementów. To zawężenie zakresu mogło wpłynąć na wyniki, szczególnie w kontekście badania skalowalności algorytmów dla naprawdę dużych instancji problemu.

Kolejnym istotnym ograniczeniem było utrzymanie stałej długości sondy (k=10) we wszystkich testach. Choć możnabyło skusić się na mniejszą wartość, tak dla większych *n* generator danych testowych zawodził, wobec tego zdecydowano o ustąpieniu w kwestii modyfikacji tego parametru. W rzeczywistych warunkach biologicznych długość ta może się znacząco różnić w zależności od zastosowanej technologii sekwencjonowania. Taka analiza mogłaby szczególnie uwidocznić różnice w zachowaniu algorytmu dokładnego i heurystycznego, gdyż im krótsze sondy, tym więcej kombinacji do rozważenia na każdym etapie rekonstrukcji.

## 5.2. Dane testowe i wyniki

W przeprowadzonych testach przyjęto następujące stałe parametry:
- Długość k-merów (sond): **k = 10**
- Stosunek błędów pozytywnych: **sqpe = n/4**

  Testy wykonane zostały wielokrotnie i uśrednione dla różnych długości sekwencji docelowych, mierząc czas wykonania obu algorytmów. Poniższa tabela przedstawia uzyskane wyniki:

| Długość sekwencji (n) | Algorytm dokładny - Czas [s] | Algorytm heurystyczny - Czas [s] | Przykładowa sekwencja (fragment)       |
|----------------------|----------------------------|---------------------------------|----------------------------------------|
| 20                   | 0.000004                   | 0.000012                        | ATGTTGTAAACCAAGATCTA                   |
| 30                   | 0.000026                   | 0.000045                        | GGGCCTTGTCTTTGCCTACTGAAGATATAT         |
| 50                   | 0.000063                   | 0.000091                        | GGGAGTACAGGCTTGAGCATATTAAAAAGACTG...   |
| 100                  | 0.000196                   | 0.000208                        | CCATCTTGCTCTCTCGATGCTGCCAGGCGTAGA...   |
| 150                  | 0.000279                   | 0.000253                        | CCCGTCTCTACTAAAAATACAAAATCAGCCAAG...   |
| 250                  | 0.000481                   | 0.000327                        | ATCTGGAGGGAAAACTGTGTGGAGAGAACACT...   |
| 500                  | 0.001098                   | 0.000782                        | AAGTGGGAGAAAAAGCTGCTGCCCATCCAGCA...   |
| 1000                 | 0.002375                   | 0.001463                        | AATCAACTCAAAATTGATTGAACATGTAAATA...   |   

Dla bardzo krótkich sekwencji algorytm dokładny wykazuje nieznaczną przewagę czasową.
 

## 5.3. Wnioski

W przeprowadzonych eksperymentach zaobserwowano, że w niektórych przypadkach algorytm dokładny (BFS) osiągał lepsze czasy wykonania niż algorytm heurystyczny (Beam Search), mimo teoretycznych założeń sugerujących odwrotną zależność. Zjawisko to można wyjaśnić poprzez analizę charakterystyki danych wejściowych. Im mniejszego kalibru były dane, tym korzystniej wypadał BFS. Beam Search, pomimo teoretycznej przewagi, traci czas na zarządzanie strukturą kopca i obliczanie heurystyk, co w małych przestrzeniach stanów staje się niepotrzebnym obciążeniem. Również, dla małych instancji problemu przestrzeń stanów w BFS jest na tyle ograniczona, że pełne przeszukiwanie nie generuje znaczącego narzutu obliczeniowego.

Dla większych problemów warto zastosować Beam Search, ale konieczne może być dostrojenie heurystyki - tzn. znalezienie innego czynnika wpływającego na ocenę kandydata w celu dalszej analizy. Dla małych problemów, gdzie przestrzeń stanów jest ograniczona, BFS jest wystarczający i często szybszy dzięki prostocie implementacji oraz statystycznie mniejszej złożoność pamięci w tej instancji problemu.

---

