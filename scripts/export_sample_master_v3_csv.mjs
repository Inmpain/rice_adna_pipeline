import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "/Users/inmpain/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const input = process.argv[2] ?? "/Users/inmpain/Desktop/rice/outputs/01a08606-6807-7130-afd5-dcca0e546e47/rice_environmental_sample_master_v2.xlsx";
const output = process.argv[3] ?? "/Users/inmpain/Desktop/rice_adna_pipeline/manifests/rice_environmental_sample_master_v3.csv";

// Coordinates reported by exact CAM core in sample_meta_data_20250922.tsv.
const angkorCoreCoordinates = new Map([
  ["CAM23-11", [13.412944, 103.864806]],
  ["CAM23-13", [13.413000, 103.864639]],
  ["CAM22-08", [13.427250, 103.844722]],
  ["CAM2201", [13.427556, 103.766056]],
  ["CAM2509", [13.412980, 103.864648]],
]);

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(input));
const values = workbook.worksheets.getItem("sample_master").getUsedRange().values;
const headers = values[0].map(String);
const latIndex = headers.indexOf("latitude");
const lonIndex = headers.indexOf("longitude");
const projectIndex = headers.indexOf("project");
const siteIndex = headers.indexOf("site");
const coreIndex = headers.indexOf("core_id");
if ([latIndex, lonIndex, projectIndex, siteIndex, coreIndex].some((x) => x < 0)) {
  throw new Error("Required metadata columns are missing");
}

const outputHeaders = [...headers, "coordinate_status", "coordinate_source", "coordinate_note"];
const rows = values.slice(1).map((sourceRow) => {
  const row = [...sourceRow];
  const project = String(row[projectIndex] ?? "");
  const site = String(row[siteIndex] ?? "");
  const core = String(row[coreIndex] ?? "");
  const hasCoordinates = row[latIndex] !== null && row[latIndex] !== ""
    && row[lonIndex] !== null && row[lonIndex] !== "";
  let status = "";
  let source = "";
  let note = "";

  if (project === "angkor" && angkorCoreCoordinates.has(core)) {
    const [lat, lon] = angkorCoreCoordinates.get(core);
    row[latIndex] = lat;
    row[lonIndex] = lon;
    status = "reported_in_cam_metadata";
    source = "sample_meta_data_20250922.tsv: field_sample_id@field_sample, latitude, longitude";
    note = `Exact coordinate reported for Angkor core ${core}; independently reproduced in angkor_final_metadata.tsv`;
  } else if (hasCoordinates) {
    status = "reported_in_sample_metadata";
    source = "nanzuo_mcp_sample_metadata.tsv";
    note = "Coordinates copied unchanged from source sample metadata";
  } else if (!site) {
    status = "missing_site_metadata";
    note = "Site is unknown; latitude and longitude cannot be assigned safely";
  } else {
    status = "needs_georeferencing";
    note = "Site name is available, but no verified sampling-point coordinate is present";
  }
  return [...row, status, source, note];
});

const escapeCsv = (value) => {
  if (value === null || value === undefined) return "";
  const text = String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
};
const csv = [outputHeaders, ...rows].map((row) => row.map(escapeCsv).join(",")).join("\n") + "\n";
await fs.writeFile(output, csv, "utf8");

const counts = {};
for (const row of rows) counts[row[headers.length]] = (counts[row[headers.length]] || 0) + 1;
console.log(JSON.stringify({ output, rows: rows.length, columns: outputHeaders.length, coordinateStatus: counts }, null, 2));
