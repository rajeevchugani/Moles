# Mole Lab

Interactive IGCSE mole calculation revision site. It is plain HTML, CSS and JavaScript.
There is nothing to build or install, and no database or server-side code.

## What it does

Mole Lab helps students practise IGCSE mole calculations by doing them, not just reading them.
Everything runs in the browser. Nothing is saved or sent anywhere.

### The tabs

1. **Reaction lab** (the main activity). Five steps:
   1. Pick a reaction type (synthesis, combustion, decomposition, displacement, acid reactions: metal, metal oxide, metal hydroxide or metal carbonate + acid) and a specific reaction.
   2. Read the word equation and write the formula of each substance.
   3. Balance the equation with the +/- coefficient buttons.
   4. Set the amounts and use the virtual bench: pour, weigh on a balance, collect gas in a syringe, or run a titration.
   5. Work through the calculation grid under the equation: Mr, moles, moles / coefficient, limiting and excess reagent, product amounts, concentration and % yield.
2. **Balancing**: balance equations with instant checking.
3. **Mole calculations**: random questions on mass, concentration and gas volume.
4. **% yield**: stand-alone percentage yield questions.
5. **Spot the mistake**: a worked solution with one error in it. Students click the first wrong line, choose the type of error, then type the corrected answer.
6. **Test**: ten mixed questions with no hints. Marks and full working are shown at the end, with the topics to revise.

### Extras

- **Streak**: a counter of correct answers in a row. After every five, students are offered the next difficulty mode.
- **Light / Dark switch** in the header. The choice, the difficulty mode and the last tab are remembered in each student's browser.
- **Link back** to the Science Revision Hub at the top of the page (edit the `hublink` address in `index.html` if the hub moves).

### Difficulty modes (top of the page)

- **Worked**: answers are shown and explained. Good for first exposure or a demonstration.
- **Guided**: fewer blanks, instant right/wrong feedback and hints.
- **Standard**: more blanks, check when ready, hints on request.
- **Challenge**: most blanks, no hints.

### Tools in the header

- **Calculator**: a floating calculator that can be dragged anywhere, minimised or switched off.
- **Periodic table**: shows atomic number, symbol and relative atomic mass. Hover or click an element
  (click keeps it pinned). "Select several" builds a short list with repeated elements and a total Mr,
  and a formula can be typed in to get its Mr. It can be resized, faded so students can see and type through it, or minimised.

### Suggested use in class

- Demonstrate one reaction in **Worked** mode on the board.
- Students move to **Guided**, then **Standard**, on their own devices.
- Use **Spot the mistake** as a starter or plenary; it works well in pairs.
- Use **Test** as a short plenary or revision homework.
- Set **Challenge** mode for extension or revision homework.

## Files

- `index.html`  the page
- `css/style.css`  all styling (dark navy theme matching the Science Revision Hub)
- `js/app.js`  all the chemistry logic and the interface
- `Dockerfile`  optional, serves the folder with nginx

## Run it

**Any web server.** Upload the whole `mole-lab` folder (keeping the folder structure) to your
web root or a sub-folder, for example `https://yourschool.org/mole-lab/`.
Works on Apache, nginx, IIS, shared hosting, GitHub Pages, Netlify, school intranets and so on.

**Quick test on your own computer.** In a terminal inside the `mole-lab` folder:

    python3 -m http.server 8000

then open http://localhost:8000

**Docker.**

    docker build -t mole-lab .
    docker run -p 8080:80 mole-lab

then open http://localhost:8080

Opening `index.html` directly by double-click also works for most features.

## Notes

- The header has on/off buttons for a floating calculator and a periodic table. Both can be dragged
  and minimised. The periodic table has Small/Medium/Large sizes, a Fade option and a Compact mode (search box only). It shows atomic number, symbol and relative atomic mass, with search.

- Fonts (Bricolage Grotesque, DM Sans, DM Mono) load from Google Fonts. If the
  network blocks Google, the site falls back to the device's system fonts and still works.
  To host the fonts yourself, download them and replace the Google Fonts link in `index.html`
  with `@font-face` rules in `css/style.css`.
- The calculator position and on/off setting are remembered in each student's browser
  (localStorage). Nothing is sent anywhere. The site collects no data.
- To change the reactions, edit the lists near the top of `js/app.js`
  (`DEC`, `FUELS`, `ACIDS`, `BASES`, `NBASE`, `DPART`). Equations are balanced automatically.
