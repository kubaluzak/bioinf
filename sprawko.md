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

## 2. Algorytm BFS — dokładny

Do rozwiązania problemu wykorzystujemy algorytm **przeszukiwania wszerz (BFS)** po przestrzeni stanów reprezentującej możliwe rozszerzenia rekonstrukcji sekwencji.

### 2.1. Reprezentacja stanu

Każdy stan algorytmu reprezentowany jest przez:

- **Suffix sekwencji** z ostatnich \(k-1\) nukleotydów zakodowanych w spektrum WS,  
- **Suffix sekwencji** z ostatnich \(k-1\) nukleotydów zakodowanych w spektrum RY,  
- Aktualnie zrekonstruowaną sekwencję nukleotydów (ciąg liter A, C, T, G).

W kodzie, zamiast przechowywać całe sekwencje, operujemy na suffixach i rozszerzamy je o kolejne nukleotydy zgodne z sondami.

### 2.2. Inicjalizacja BFS

- Na początku kodujemy suffix startowej sekwencji w obu spektrach,  
- Tworzymy dwie mapy prefixów (po jednym na każde spektrum), które dla prefixu długości \(k-1\) podają zbiór możliwych nukleotydów do dopisania,  
- Inicjujemy kolejkę BFS z początkowym stanem: suffixy startowe i startowa sekwencja,  
- Tworzymy zbiór `visited` by pamiętać, które pary suffixów już zostały przetworzone — to zapobiega powtórkom.

### 2.3. Iteracyjny proces BFS

W każdej iteracji:

1. Pobieramy stan z kolejki — aktualne suffixy WS i RY oraz dotychczas zrekonstruowaną sekwencję.  
2. Sprawdzamy, czy długość sekwencji osiągnęła wymaganą wartość — jeśli tak, to zwracamy rozwiązanie.  
3. Pobieramy z map prefixów kandydujące nukleotydy, które można dołączyć do aktualnych suffixów WS i RY.  
4. Iterujemy przez nukleotydy \( \{A, C, T, G\} \) i sprawdzamy, czy nukleotyd pasuje jednocześnie do obu spektrów (czy występuje jako dozwolony w obu mapach dla odpowiednich prefixów).  
5. Dla każdego pasującego nukleotydu tworzymy nowy stan:  
   - Obliczamy nowe suffixy, przesuwając istniejący suffix o 1 pozycję i dopisując kodowaną literę (w spektrum WS i RY),  
   - Doklejamy nukleotyd do aktualnej sekwencji,  
   - Jeśli para suffixów nie była wcześniej odwiedzona, dodajemy nowy stan do kolejki i do zbioru `visited`.

### 2.4. Zakończenie i wynik

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

## 4. Zależność czasu działania algorytmu od długości sekwencji
Parametry dla każdego przypadku: 

k = 10

sqpe = n/4

pose = n/2

| Długość sekwencji | Czas wykonania [s] | Przykładowa sekwencja (fragment)               |
|-------------------|--------------------|-----------------------------------------------|
| 20                | 0.000004           | ATGTTGTAAACCAAGATCTA                          |
| 30                | 0.000026           | GGGCCTTGTCTTTGCCTACTGAAGATATAT                |
| 50                | 0.000063           | GGGAGTACAGGCTTGAGCATATTAAAAAGACTGCTTTTTAATAT |
| 100               | 0.000196           | CCATCTTGCTCTCTCGATGCTGCCAGGCGTAGAAGTTGAGC... |
| 150               | 0.000279           | CCCGTCTCTACTAAAAATACAAAATCAGCCAAGCGTGGTGGC... |
| 250               | 0.000481           | ATCTGGAGGGAAAACTGTGTGGAGAGAACACTTGACAAGAAA... |
| 500               | 0.001098           | AAGTGGGAGAAAAAGCTGCTGCCCATCCAGCAATGGAGCTTC... |
| 1000              | 0.002375           | AATCAACTCAAAATTGATTGAACATGTAAATATAAGACCTGA... |

---

# Heurystyczny algorytm sekwencjonowania DNA

## Opis algorytmu

Algorytm heurystyczny służy do **rekonstrukcji sekwencji DNA** z dwóch spektrów binarnych (WS i RY) generowanych metodą hybrydyzacji. Zamiast przeszukiwać wszystkie możliwe kombinacje, algorytm stosuje strategię zachłanną (greedy), która w każdej iteracji wybiera najbardziej obiecujący krok na podstawie zgodności z sondami.

---

## Założenia

- Znana jest **początkowa sekwencja startowa** o długości co najmniej \(k-1\),
- Znana jest **docelowa długość sekwencji** \(n\),
- Znane są **binarne spektra** WS i RY, które mogą zawierać błędy pozytywne (fałszywe jedynki),
- Algorytm wybiera nukleotydy pasujące do **obu spektrów jednocześnie**.

---

## Reprezentacja danych

- **Spektrum binarne** — dla każdego możliwego k-meru (ciągu długości \(k\)) istnieje informacja o jego obecności (1) lub braku (0) w spektrum WS i RY.
- **Prefix mapy** — dla każdego możliwego prefixu długości \(k-1\) przechowujemy listę możliwych dopisań na podstawie spektrum.

---

## Działanie algorytmu

### Krok 1: Inicjalizacja

1. Zakoduj początkowy **suffix** (ostatnie \(k-1\) znaków) startowej sekwencji w obu spektrach.
2. Zbuduj **mapy prefixów** WS i RY, które dla każdego prefixu podają zbiór możliwych kolejnych nukleotydów.

---

### Krok 2: Iteracyjne rozszerzanie sekwencji

Dla każdej pozycji od długości startowej do docelowej długości \(n\):

1. Pobierz listę możliwych nukleotydów do dopisania dla aktualnego suffixu w **obu** spektrach.
2. Wyznacz **część wspólną** obu zbiorów kandydatów — tylko nukleotydy akceptowalne przez oba spektra.
3. Jeśli istnieje tylko **jeden kandydat**, dopisz go do sekwencji.
4. Jeśli istnieje **więcej niż jeden kandydat**, zastosuj heurystykę wyboru:
   - Preferuj nukleotydy zgodne z kolejnymi możliwymi sondami,
   - Możesz zastosować prostą strategię np. **alfabetyczne pierwszeństwo** lub **liczbę możliwych rozszerzeń w kolejnych krokach**.
5. Jeśli **nie ma żadnego możliwego dopisania**, zakończ rekonstrukcję z wynikiem negatywnym.
6. Po dopisaniu nukleotydu zaktualizuj aktualny suffix i przejdź do kolejnego kroku.

---

### Krok 3: Wynik

- **Sukces**: zrekonstruowano sekwencję o długości \(n\),
- **Niepowodzenie**: nie znaleziono dopasowania w jednym z kroków.

---

## Pseudokod

```plaintext
HEURYSTYKA(start_sequence, spektrum_WS, spektrum_RY, k, n):
    suffix_WS ← koduj(start_sequence, spektrum_WS)
    suffix_RY ← koduj(start_sequence, spektrum_RY)
    map_WS ← zbuduj_mapę_prefixów(spektrum_WS)
    map_RY ← zbuduj_mapę_prefixów(spektrum_RY)

    sequence ← start_sequence

    while length(sequence) < n:
        candidates_WS ← map_WS[suffix_WS]
        candidates_RY ← map_RY[suffix_RY]
        candidates ← część_wspólna(candidates_WS, candidates_RY)

        if candidates is empty:
            return "Brak rozwiązania"

        next_nucleotide ← wybierz_kandydata(candidates)
        sequence ← sequence + next_nucleotide

        suffix_WS ← przesun(suffix_WS, next_nucleotide)
        suffix_RY ← przesun(suffix_RY, next_nucleotide)

    return sequence

