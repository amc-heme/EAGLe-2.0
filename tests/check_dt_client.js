// Execute DT's installed browser binding at the exact hidden-widget failure point.
const fs = require('fs');
const vm = require('vm');
const assert = require('assert');
let binding;
const context = {
  window: {},
  HTMLWidgets: {shinyMode: false, widget: value => { binding = value; }}
};
vm.runInNewContext(fs.readFileSync(process.argv[2], 'utf8'), context);
assert.strictEqual(binding.name, 'datatables');
const hidden = {offsetWidth: 0, offsetHeight: 0};
// Establish that this installed DT binding reproduces the user's exact error.
assert.throws(() => binding.renderValue(hidden, null, {}), /lazyRender/);
const events = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
for (const event of events) {
  if (event.suspended) continue; // Shiny sends no widget payload for a suspended output.
  assert.notStrictEqual(event.payload.x, null, `${event.output} emitted x:null`);
  const instance = {};
  binding.renderValue(hidden, event.payload.x, instance);
  assert.strictEqual(instance.data, event.payload.x);
}
console.log(`PASS: installed DT JavaScript reproduces the old error; ${events.length} actual output transitions avoid it`);
